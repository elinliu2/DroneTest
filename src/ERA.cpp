#include "DroneTrajectory.h"
#include "Splotting.h"

zkpk DroneTrajectory::theGigaAlgo(SystemState currState)
{
    Eigen::Vector<double, NUM_STATES> z0;
    z0 << currState.plant, currState.alge;
    Eigen::Vector<double, NUM_PARAMETERS> pk_prev = Eigen::Vector<double, NUM_PARAMETERS>::Zero();
    zkpk curr = {z0, getParams()};
    int count = 0;
    while((curr.pk-pk_prev).cwiseAbs().sum() > 1e-8)
    {
        pk_prev = curr.pk;
        curr = updateStep(curr, z0);
    }
    m_logger << "m_finalTime: " << m_finalTime << std::endl;
    return curr;
}

zkpk DroneTrajectory::updateStep(zkpk prev, Eigen::Vector<double, NUM_STATES> currState)
{
    SystemState prev_zk_state = {prev.zk.segment(0, NUM_PLANT_STATES), prev.zk.segment(NUM_PLANT_STATES, NUM_ALGE_STATES)};
    SimResults traj = Trajectory(prev_zk_state);
    std::vector<dwdwo> ts = trajSens(traj);
    G_tp gtp = calc_G_tp(ts);
    Eigen::Vector<double, NUM_STATES+NUM_PARAMETERS> dG = calc_dG_test(prev_zk_state, ts.at(gtp.tp), gtp, traj.time.at(gtp.tp) + m_simTimestep);
    // for(int i = 0; i < NUM_STATES+NUM_PARAMETERS; i++) {m_logger << dG(i) << std::endl;} m_logger << std::endl;
    Eigen::Vector<double, NUM_STATES> vz = dG.segment(0, NUM_STATES);
    Eigen::Vector<double, NUM_PARAMETERS> vp = dG.segment(NUM_STATES, NUM_PARAMETERS);
    Eigen::Vector<double, NUM_PARAMETERS> pk = prev.pk + m_algo_alpha*(vz.transpose()*(prev.zk - currState) + gtp.G - m_epsilon)/(vz.transpose()*m_Pinv*vz)*vp;
    setParams(pk);
    Eigen::Vector<double, NUM_STATES> zk = (m_Pinv*vz*vz.transpose())/(vz.transpose()*m_Pinv*vz)*(prev.zk-currState) + (m_Pinv*vz*vp.transpose())/(vz.transpose()*m_Pinv*vz)*(prev.pk-pk) + currState - (m_Pinv*vz)/(vz.transpose()*m_Pinv*vz)*(gtp.G-m_epsilon);
    double backtrack = 0.5;
    traj = Trajectory({zk.segment(0, NUM_PLANT_STATES), zk.segment(NUM_PLANT_STATES, NUM_ALGE_STATES)});
    while((!traj.stable || !traj.converged) && (pk-prev.pk).cwiseAbs().sum() > 1e-8)
    {
        pk = prev.pk + backtrack*(pk-prev.pk);
        setParams(pk);
        zk = prev.zk + backtrack*(zk-prev.zk);
        traj = Trajectory({zk.segment(0, NUM_PLANT_STATES), zk.segment(NUM_PLANT_STATES, NUM_ALGE_STATES)});
        backtrack /= 2;
    }
    if(!traj.stable || !traj.converged){
        // SimResults traj = Trajectory(prev_zk_state);
        // m_logger << "traj converge?: " << traj.converged << " stable?: " << traj.stable << std::endl;
        // SystemState z0_state = {currState.segment(0, NUM_PLANT_STATES), currState.segment(NUM_PLANT_STATES, NUM_ALGE_STATES)};
        // SimResults z0_state_traj = Trajectory(z0_state);
        // m_logger << "z0_state traj converge?: " << z0_state_traj.converged << " stable?: " << z0_state_traj.stable << std::endl;
        return prev;
    } else {
        // m_logger << "traj converge?: " << traj.converged << " stable?: " << traj.stable << std::endl;
        // SystemState z0_state = {currState.segment(0, NUM_PLANT_STATES), currState.segment(NUM_PLANT_STATES, NUM_ALGE_STATES)};
        // SimResults z0_state_traj = Trajectory(z0_state);
        // m_logger << "z0_state traj converge?: " << z0_state_traj.converged << " stable?: " << z0_state_traj.stable << std::endl;
        return {zk, pk};
    }
}

Eigen::Vector<double, NUM_STATES> DroneTrajectory::closestZBar(SystemState currState)
{
    Eigen::Vector<double, NUM_STATES> z0;
    z0 << currState.plant, currState.alge;

    Eigen::Vector<double, NUM_STATES> zk = z0;
    Eigen::Vector<double, NUM_STATES> zk_prev = z0;
    do {
        SimResults traj = Trajectory({zk.segment(0, NUM_PLANT_STATES), zk.segment(NUM_PLANT_STATES, NUM_ALGE_STATES)});
        std::vector<dwdwo> ts = trajSens(traj);
        G_tp gtp = calc_G_tp(ts);
        std::vector<dwdp> ts_p = trajSensParam(traj, gtp);
        Eigen::Vector<double, NUM_STATES+NUM_PARAMETERS> dG = calc_dG_test({zk.segment(0, NUM_PLANT_STATES), zk.segment(NUM_PLANT_STATES, NUM_ALGE_STATES)}, ts.at(gtp.tp), gtp, traj.time.at(gtp.tp));
        Eigen::Vector<double, NUM_STATES> vz = dG.segment(0, NUM_STATES);

        zk_prev = zk;
        zk = z0 + (m_Pinv*vz*vz.transpose())/(vz.transpose()*m_Pinv*vz)*(zk - z0) - (m_Pinv*vz)/(vz.transpose()*m_Pinv*vz)*(gtp.G - m_epsilon);

        double backtrack = 0.1;
        traj = Trajectory({zk.segment(0, NUM_PLANT_STATES), zk.segment(NUM_PLANT_STATES, NUM_ALGE_STATES)});
        while(!traj.stable && (zk-zk_prev).cwiseAbs().sum() > 1e-8)
        {
            zk = zk_prev + backtrack*(zk-zk_prev);
            traj = Trajectory({zk.segment(0, NUM_PLANT_STATES), zk.segment(NUM_PLANT_STATES, NUM_ALGE_STATES)});
            backtrack /= 2;
        }
        if(!traj.stable){
            zk = zk_prev;
        }
    } while((zk-zk_prev).cwiseAbs().sum() > 1e-8);

    return zk;
}

Eigen::Vector<double, NUM_PARAMETERS> DroneTrajectory::getParams()
{
    Eigen::Vector<double, NUM_PARAMETERS> params;

    params(kppx) = m_ctrlParams.at(posX).kp;
    params(kipx) = m_ctrlParams.at(posX).ki;
    params(kdpx) = m_ctrlParams.at(posX).kd;
    params(kppy) = m_ctrlParams.at(posY).kp;
    params(kipy) = m_ctrlParams.at(posY).ki;
    params(kdpy) = m_ctrlParams.at(posY).kd;
    params(kppz) = m_ctrlParams.at(posZ).kp;
    params(kipz) = m_ctrlParams.at(posZ).ki;
    params(kdpz) = m_ctrlParams.at(posZ).kd;

    params(kpvx) = m_ctrlParams.at(velX).kp;
    params(kivx) = m_ctrlParams.at(velX).ki;
    params(kdvx) = m_ctrlParams.at(velX).kd;
    params(kpvy) = m_ctrlParams.at(velY).kp;
    params(kivy) = m_ctrlParams.at(velY).ki;
    params(kdvy) = m_ctrlParams.at(velY).kd;
    params(kpvz) = m_ctrlParams.at(velZ).kp;
    params(kivz) = m_ctrlParams.at(velZ).ki;
    params(kdvz) = m_ctrlParams.at(velZ).kd;

    params(kpphi) = m_ctrlParams.at(roll).kp;
    params(kiphi) = m_ctrlParams.at(roll).ki;
    params(kdphi) = m_ctrlParams.at(roll).kd;
    params(kptheta) = m_ctrlParams.at(pitch).kp;
    params(kitheta) = m_ctrlParams.at(pitch).ki;
    params(kdtheta) = m_ctrlParams.at(pitch).kd;
    params(kppsi) = m_ctrlParams.at(yaw).kp;
    params(kipsi) = m_ctrlParams.at(yaw).ki;
    params(kdpsi) = m_ctrlParams.at(yaw).kd;

    params(kpp) = m_ctrlParams.at(rollRate).kp;
    params(kip) = m_ctrlParams.at(rollRate).ki;
    params(kdp) = m_ctrlParams.at(rollRate).kd;
    params(kpq) = m_ctrlParams.at(pitchRate).kp;
    params(kiq) = m_ctrlParams.at(pitchRate).ki;
    params(kdq) = m_ctrlParams.at(pitchRate).kd;
    params(kpr) = m_ctrlParams.at(yawRate).kp;
    params(kir) = m_ctrlParams.at(yawRate).ki;
    params(kdr) = m_ctrlParams.at(yawRate).kd;

    return params;
}

void DroneTrajectory::setParams(Eigen::Vector<double, NUM_PARAMETERS> params)
{
    m_ctrlParams.at(posX).kp = params(kppx);
    m_ctrlParams.at(posX).ki = params(kipx);
    m_ctrlParams.at(posX).kd = params(kdpx);
    m_ctrlParams.at(posY).kp = params(kppy);
    m_ctrlParams.at(posY).ki = params(kipy);
    m_ctrlParams.at(posY).kd = params(kdpy);
    m_ctrlParams.at(posZ).kp = params(kppz);
    m_ctrlParams.at(posZ).ki = params(kipz);
    m_ctrlParams.at(posZ).kd = params(kdpz);

    m_ctrlParams.at(velX).kp = params(kpvx);
    m_ctrlParams.at(velX).ki = params(kivx);
    m_ctrlParams.at(velX).kd = params(kdvx);
    m_ctrlParams.at(velY).kp = params(kpvy);
    m_ctrlParams.at(velY).ki = params(kivy);
    m_ctrlParams.at(velY).kd = params(kdvy);
    m_ctrlParams.at(velZ).kp = params(kpvz);
    m_ctrlParams.at(velZ).ki = params(kivz);
    m_ctrlParams.at(velZ).kd = params(kdvz);

    m_ctrlParams.at(roll).kp = params(kpphi);
    m_ctrlParams.at(roll).ki = params(kiphi);
    m_ctrlParams.at(roll).kd = params(kdphi);
    m_ctrlParams.at(pitch).kp = params(kptheta);
    m_ctrlParams.at(pitch).ki = params(kitheta);
    m_ctrlParams.at(pitch).kd = params(kdtheta);
    m_ctrlParams.at(yaw).kp = params(kppsi);
    m_ctrlParams.at(yaw).ki = params(kipsi);
    m_ctrlParams.at(yaw).kd = params(kdpsi);

    m_ctrlParams.at(rollRate).kp = params(kpp);
    m_ctrlParams.at(rollRate).ki = params(kip);
    m_ctrlParams.at(rollRate).kd = params(kdp);
    m_ctrlParams.at(pitchRate).kp = params(kpq);
    m_ctrlParams.at(pitchRate).ki = params(kiq);
    m_ctrlParams.at(pitchRate).kd = params(kdq);
    m_ctrlParams.at(yawRate).kp = params(kpr);
    m_ctrlParams.at(yawRate).ki = params(kir);
    m_ctrlParams.at(yawRate).kd = params(kdr);
}