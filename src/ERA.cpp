#include "DroneTrajectory.h"
#include "Splotting.h"
#include <iomanip>

zkpk DroneTrajectory::theGigaAlgo(SystemState currState)
{
    Eigen::Vector<double, NUM_STATES> z0;
    z0 << currState.plant, currState.alge;
    Eigen::Vector<double, NUM_PARAMETERS> pk_prev = Eigen::Vector<double, NUM_PARAMETERS>::Zero();
    zkpk curr = {z0, m_sf};
    Eigen::Vector<double, NUM_PARAMETERS> p0 = m_sf;
    int count = 0;
    while((curr.pk-pk_prev).cwiseAbs().sum() > 5e-5)
    {
        pk_prev = curr.pk;
        curr = updateStep(curr, z0);
        SimResults traj = Trajectory(currState);
        std::vector<dwdwo> ts = trajSens(traj);
        G_tp gtp = calc_G_tp(ts);
        count++;
        m_logger << "count " << count << std::endl;
        m_logger << "(curr.pk-pk_prev).cwiseAbs().sum() " << (curr.pk-pk_prev).cwiseAbs().sum() << std::endl;
        m_logger << "(curr.pk-p0).cwiseAbs().sum() " <<  (curr.pk-p0).cwiseAbs().sum() << std::endl;
        m_logger << "stable?: " << traj.stable  << " converged?: " << traj.converged << " G: " << gtp.G << std::endl;
    }
    m_logger << "count: " << count << std::endl;
    return curr;
}

zkpk DroneTrajectory::updateStep(zkpk prev, Eigen::Vector<double, NUM_STATES> const & currState)
{
    SystemState prev_zk_state = {prev.zk.segment(0, NUM_PLANT_STATES), prev.zk.segment(NUM_PLANT_STATES, NUM_ALGE_STATES)};
    SimResults traj = Trajectory(prev_zk_state);
    std::vector<dwdwo> ts = trajSens(traj);
    G_tp gtp = calc_G_tp(ts);
    m_logger << "tp: " << gtp.tp << std::endl;
    d2w d2w = calc_d2w(traj, ts, gtp);
    Eigen::Vector<double, NUM_STATES+NUM_PARAMETERS> dG = calc_dG(ts.at(gtp.tp), d2w, gtp);
    Eigen::Vector<double, NUM_STATES> vz = dG.segment(0, NUM_STATES);
    m_logger << "vz " << vz << std::endl;
    Eigen::Vector<double, NUM_PARAMETERS> vp = dG.segment(NUM_STATES, NUM_PARAMETERS);
    m_logger << "vp " << vp << std::endl;
    double denominator = vz.transpose()*m_Pinv*vz;
    Eigen::Vector<double, NUM_PARAMETERS> pk = prev.pk + m_algo_alpha*(vz.transpose()*(currState - prev.zk) + gtp.G - m_epsilon)/denominator*vp;
    m_sf = pk;
    Eigen::Vector<double, NUM_STATES> zk = (m_Pinv*vz*vz.transpose())/denominator*(prev.zk-currState) - (m_Pinv*vz*vp.transpose())/denominator*(pk - prev.pk) + currState - (m_Pinv*vz)/denominator*(gtp.G-m_epsilon);
    double backtrack = m_backtrack;
    traj = Trajectory({zk.segment(0, NUM_PLANT_STATES), zk.segment(NUM_PLANT_STATES, NUM_ALGE_STATES)});
    int backtrackingCount = 0;
    std::chrono::time_point start = std::chrono::steady_clock::now();
    while((!traj.stable || !traj.converged) && (pk-prev.pk).cwiseAbs().sum() > 5e-5)
    {
        m_logger << "backtrackingCount " << backtrackingCount << std::endl;
        pk = prev.pk + backtrack*(pk-prev.pk);
        m_sf = pk; 
        zk = prev.zk + backtrack*(zk-prev.zk);
        m_logger << "pk " << pk << std::endl;
        traj = Trajectory({zk.segment(0, NUM_PLANT_STATES), zk.segment(NUM_PLANT_STATES, NUM_ALGE_STATES)});
        backtrack /= 2;
        backtrackingCount++;
    }
    if(backtrackingCount > 0)
    {
        std::chrono::time_point end = std::chrono::steady_clock::now();
        std::chrono::microseconds elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
        m_logger << "backtracking time: " << elapsed.count() << " us" << std::endl;
        m_logger << "backtracking count: " << backtrackingCount << std::endl;
        if(backtrackingCount > 2)
        {
            m_backtrack /= 2.0;
        }
    }
    if(!traj.stable || !traj.converged){
        return prev;
    } else {
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
        d2wdwo2 d2wdwo2 = calc_d2wdwo2(traj, ts, gtp);
        Eigen::Vector<double, NUM_STATES> vz = calc_vz(ts.at(gtp.tp), d2wdwo2, gtp);

        zk_prev = zk;
        zk = z0 + (m_Pinv*vz*vz.transpose())/(vz.transpose()*m_Pinv*vz)*(zk - z0) - (m_Pinv*vz)/(vz.transpose()*m_Pinv*vz)*(gtp.G - m_epsilon);

        double backtrack = 0.1;
        traj = Trajectory({zk.segment(0, NUM_PLANT_STATES), zk.segment(NUM_PLANT_STATES, NUM_ALGE_STATES)});
        while(!traj.stable && (zk-zk_prev).cwiseAbs().sum() > 5e-5)
        {
            m_logger << "traj.stable " << traj.stable << " (zk-zk_prev).cwiseAbs().sum(): " << (zk-zk_prev).cwiseAbs().sum() << std::endl;
            zk = zk_prev + backtrack*(zk-zk_prev);
            traj = Trajectory({zk.segment(0, NUM_PLANT_STATES), zk.segment(NUM_PLANT_STATES, NUM_ALGE_STATES)});
            backtrack /= 2;
        }
        if(!traj.stable){
            zk = zk_prev;
        }
        m_logger << "traj.stable " << traj.stable << " (zk-zk_prev).cwiseAbs().sum(): " << (zk-zk_prev).cwiseAbs().sum() << std::endl;
    } while((zk-zk_prev).cwiseAbs().sum() > 1e-5);

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