#include "DroneTrajectory.h"
#include "Splotting.h"
#include <iomanip>

zkpk DroneTrajectory::theGigaAlgo(SystemState currState)
{
    Eigen::Vector<double, NUM_STATES> z0;
    z0 << currState.plant, currState.alge;
    Eigen::Vector<double, NUM_PARAMETERS> pk_prev = Eigen::Vector<double, NUM_PARAMETERS>::Zero();
    zkpk curr = {z0, getParams()};
    Eigen::Vector<double, NUM_PARAMETERS> p0 = getParams();
    int count = 0;
    while((curr.pk-pk_prev).cwiseAbs().sum() > 1e-8)
    {
        pk_prev = curr.pk;
        curr = updateStep(curr, z0);
        count++;
        m_logger << "count " << count << std::endl;
        m_logger << "(curr.pk-pk_prev).cwiseAbs().sum() " << (curr.pk-pk_prev).cwiseAbs().sum() << std::endl;
        m_logger << "(curr.pk-p0).cwiseAbs().sum() " <<  (curr.pk-p0).cwiseAbs().sum() << std::endl;
    }
    m_logger << "count: " << count << std::endl;
    return curr;
}

zkpk DroneTrajectory::theGigaAlgopt2(SystemState currState)
{
    Eigen::Vector<double, NUM_STATES> z0;
    z0 << currState.plant, currState.alge;
    Eigen::Vector<double, NUM_PARAMETERS> pk_prev = Eigen::Vector<double, NUM_PARAMETERS>::Zero();
    zkpk curr = {z0, getParams()};
    Eigen::Vector<double, NUM_PARAMETERS> p0 = getParams();
    int count = 0;
    while((curr.pk-pk_prev).cwiseAbs().sum() > 1e-8)
    {
        pk_prev = curr.pk;
        curr = updateStepWNewBackstepping(curr, z0);
        count++;
        m_logger << "count " << count << std::endl;
        m_logger << "(curr.pk-pk_prev).cwiseAbs().sum() " << (curr.pk-pk_prev).cwiseAbs().sum() << std::endl;
        m_logger << "(curr.pk-p0).cwiseAbs().sum() " <<  (curr.pk-p0).cwiseAbs().sum() << std::endl;
    }
    m_logger << "total count: " << count << std::endl;
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
    setParams(pk);
    Eigen::Vector<double, NUM_STATES> zk = (m_Pinv*vz*vz.transpose())/denominator*(prev.zk-currState) - (m_Pinv*vz*vp.transpose())/denominator*(pk - prev.pk) + currState - (m_Pinv*vz)/denominator*(gtp.G-m_epsilon);
    double backtrack = m_backtrack;
    traj = Trajectory({zk.segment(0, NUM_PLANT_STATES), zk.segment(NUM_PLANT_STATES, NUM_ALGE_STATES)});
    int backtrackingCount = 0;
    std::chrono::time_point start = std::chrono::steady_clock::now();
    while((!traj.stable || !traj.converged) && (pk-prev.pk).cwiseAbs().sum() > 1e-8)
    {
        m_logger << "backtrackingCount " << backtrackingCount << std::endl;
        pk = prev.pk + backtrack*(pk-prev.pk);
        setParams(pk);
        zk = prev.zk + backtrack*(zk-prev.zk);
        m_logger << "pk " << pk << std::endl;
        traj = Trajectory({zk.segment(0, NUM_PLANT_STATES), zk.segment(NUM_PLANT_STATES, NUM_ALGE_STATES)});
        backtrack /= 2;
        backtrackingCount++;
        // m_logger << "sim time length: " << traj.time.size() << std::endl;
        // if (traj.stable && traj.converged)
        // {
        //     m_logger << "final state: " << std::endl;
        //     m_logger << "x: " << traj.stateProgression.at(traj.time.size()-1).plant(x) <<
        //                 " y: " << traj.stateProgression.at(traj.time.size()-1).plant(y) <<
        //                 " z: " << traj.stateProgression.at(traj.time.size()-1).plant(z) <<
        //                 " phi: " << traj.stateProgression.at(traj.time.size()-1).plant(phi) <<
        //                 " theta: " << traj.stateProgression.at(traj.time.size()-1).plant(theta) <<
        //                 " psi: " << traj.stateProgression.at(traj.time.size()-1).plant(psi) <<
        //                 " xdot: " << traj.stateProgression.at(traj.time.size()-1).plant(xdot) <<
        //                 " ydot: " << traj.stateProgression.at(traj.time.size()-1).plant(ydot) <<
        //                 " zdot: " << traj.stateProgression.at(traj.time.size()-1).plant(zdot) <<
        //                 " p: " << traj.stateProgression.at(traj.time.size()-1).plant(p) <<
        //                 " q: " << traj.stateProgression.at(traj.time.size()-1).plant(q) <<
        //                 " r: " << traj.stateProgression.at(traj.time.size()-1).plant(r) << std::endl;
        // }
    }
    if(backtrackingCount > 0)
    {
        std::chrono::time_point end = std::chrono::steady_clock::now();
        std::chrono::microseconds elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
        m_logger << "backtracking time: " << elapsed.count() << " us" << std::endl;
        m_logger << "backtracking count: " << backtrackingCount << std::endl;
        if(backtrackingCount > 2)
        {
            m_backtrack = backtrack;
        }
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

Eigen::Vector<double, NUM_STATES> DroneTrajectory::backstep_btwn_zk_zkp1(Eigen::Vector<double, NUM_STATES> const & zk, Eigen::Vector<double, NUM_STATES> const & zkp1)
{
    double backstep = 0.5;
    Eigen::Vector<double, NUM_STATES> intermediate_zk = backstep*zkp1 + (1-backstep)*zk;
    SystemState intermediate_zk_state = {intermediate_zk.segment(0, NUM_PLANT_STATES), intermediate_zk.segment(NUM_PLANT_STATES, NUM_ALGE_STATES)};
    SimResults traj = Trajectory(intermediate_zk_state);
    while(!traj.stable || !traj.converged){
        backstep /= 2;
        intermediate_zk = backstep*zkp1 + (1-backstep)*zk;
        intermediate_zk_state = {intermediate_zk.segment(0, NUM_PLANT_STATES), intermediate_zk.segment(NUM_PLANT_STATES, NUM_ALGE_STATES)};
        traj = Trajectory(intermediate_zk_state);
    }
    return intermediate_zk;
}

Eigen::Vector<double, NUM_STATES> vz_stateCloseToRoABoundary_testERA()
{
    Eigen::Vector<double, NUM_STATES> vz = {
        -2.84488e-08,
         6.49673e-11,
        -1.27123e-08,
        -3.72008e-09,
        -7.05403e-08,
        -1.06216e-07,
        -1.36627e-07,
         6.36192e-10,
        -5.25603e-08,
        -6.10231e-11,
         1.81503e-10,
        -8.21654e-09,
                  -0,
                  -0,
                  -0,
                  -0,
                  -0,
                  -0,
           1.285e-09,
                  -0,
                  -0,
         4.73351e-10,
                  -0,
        -6.61786e-14,
                  -0,
        -4.24746e-11,
                  -0,
                  -0,
        -7.76228e-11,
                  -0,
                  -0,
        -1.63984e-08,
                  -0,
                  -0,
        -3.92802e-10,
                  -0,
                  -0,
         6.09986e-16,
        -3.16742e-16,
         3.98727e-15,
        -2.25276e-15,
        -3.78777e-11,
                  -0,
                  -0,
        -1.74299e-09,
                  -0,
                  -0,
        -5.48197e-11,
                  -0,
                  -0,
                  -0,
                  -0,
                  -0,
                  -0,
        -4.30932e-11,
        -5.24444e-09,
        -2.58058e-10,
        -2.79594e-07,
                  -0,
                  -0
    };

    return vz;
}

Eigen::Vector<double, NUM_PARAMETERS> vp_stateCloseToRoABoundary_testERA()
{
    Eigen::Vector<double, NUM_PARAMETERS> vp = {
         -3.4627e-07,
        -4.68415e-06,
         5.17365e-07,
        -3.50998e-11,
         6.89135e-12,
        -1.65925e-11,
          1.7199e-10,
         1.35183e-11,
         1.13221e-09,
         7.49887e-09,
        -3.66293e-07,
         8.67092e-09,
        -2.87432e-12,
        -3.69854e-13,
         2.07048e-11,
         5.56554e-11,
         5.95553e-12,
         3.12149e-10,
         -9.8693e-11,
        -7.71553e-11,
         9.20922e-10,
        -5.45863e-09,
         2.92891e-08,
        -5.88542e-08,
         5.27344e-10,
         5.03168e-10,
        -5.56355e-10,
         2.27796e-13,
         4.28667e-15,
        -1.75006e-11,
        -1.60923e-11,
        -4.95524e-13,
         3.39489e-10,
         2.46575e-11,
         7.73515e-12,
        -2.13172e-11
    };
    return vp;
}


zkpk DroneTrajectory::updateStepWNewBackstepping(zkpk prev, Eigen::Vector<double, NUM_STATES> const & currState)
{
    SystemState prev_zk_state = {prev.zk.segment(0, NUM_PLANT_STATES), prev.zk.segment(NUM_PLANT_STATES, NUM_ALGE_STATES)};
    SimResults traj = Trajectory(prev_zk_state);
    std::vector<dwdwo> ts = trajSens(traj);
    G_tp gtp = calc_G_tp(ts);
    m_logger << "G: " << gtp.G <<  " tp: " << gtp.tp << std::endl;
    d2w d2w = calc_d2w(traj, ts, gtp);
    Eigen::Vector<double, NUM_STATES+NUM_PARAMETERS> dG = calc_dG(ts.at(gtp.tp), d2w, gtp);
    Eigen::Vector<double, NUM_STATES> vz = dG.segment(0, NUM_STATES);
    // Eigen::Vector<double, NUM_STATES> vz = vz_stateCloseToRoABoundary_testERA();
    m_logger << std::setprecision(15) << "vz " << vz << std::endl;
    Eigen::Vector<double, NUM_PARAMETERS> vp = dG.segment(NUM_STATES, NUM_PARAMETERS);
    // Eigen::Vector<double, NUM_PARAMETERS> vp = vp_stateCloseToRoABoundary_testERA();
    m_logger << std::setprecision(15) << "vp " << vp << std::endl;
    double denominator = vz.transpose()*m_Pinv*vz;
    Eigen::Vector<double, NUM_PARAMETERS> pk = prev.pk + m_algo_alpha*(vz.transpose()*(currState - prev.zk) + gtp.G - m_epsilon)/denominator*vp;
    // Eigen::Vector<double, NUM_PARAMETERS> pk = prev.pk + m_algo_alpha*vp;
    m_logger << "m_algo_alpha*(vz.transpose()*(currState - prev.zk) + gtp.G - m_epsilon) " << m_algo_alpha*(vz.transpose()*(currState - prev.zk) + gtp.G - m_epsilon) << std::endl;
    m_logger << "denominator " << denominator << std::endl;
    m_logger << "pk " << pk << std::endl;
    setParams(pk);
    Eigen::Vector<double, NUM_STATES> zk = (m_Pinv*vz*vz.transpose())/denominator*(prev.zk-currState) - (m_Pinv*vz*vp.transpose())/denominator*(pk - prev.pk) + currState - (m_Pinv*vz)/denominator*(gtp.G-m_epsilon);
    traj = Trajectory({zk.segment(0, NUM_PLANT_STATES), zk.segment(NUM_PLANT_STATES, NUM_ALGE_STATES)});
    int backtrackingCount = 0;
    std::chrono::time_point start = std::chrono::steady_clock::now();
    double alpha_backstep = m_algo_alpha/2;
    if (traj.stable && traj.converged) {
        m_logger << "traj stable and converged" << std::endl;
    }
    while((!traj.stable || !traj.converged) && (pk-prev.pk).cwiseAbs().sum() > 1e-8)
    {
        m_logger << "traj not stable or not converged" << std::endl;
        SimResults prev_zk_pkp1_traj = Trajectory(prev_zk_state);
        if(prev_zk_pkp1_traj.stable){
            m_logger << "prev traj stable" << std::endl;
            std::vector<dwdwo> prev_zk_pkp1_ts = trajSens(prev_zk_pkp1_traj);
            G_tp prev_zk_pkp1_gtp = calc_G_tp(prev_zk_pkp1_ts);
            m_logger << std::setprecision(15) << "prev_zk_pkp1_gtp.G " << prev_zk_pkp1_gtp.G << " G " << gtp.G << std::endl;  
            if(prev_zk_pkp1_gtp.G > gtp.G){
                // zk = closestZBar(prev_zk_state);
                m_logger << "backstep btwn zk zkp1" << std::endl;
                zk = backstep_btwn_zk_zkp1(prev.zk, zk);
                traj = Trajectory({zk.segment(0, NUM_PLANT_STATES), zk.segment(NUM_PLANT_STATES, NUM_ALGE_STATES)});
                break;
            }
            else {
                m_logger << "backstepping parameters stable" << std::endl;
                pk = prev.pk + alpha_backstep*vp;
                setParams(pk);
                alpha_backstep/=2.0;    
            }
        } else {
            m_logger << "backstepping parameters unstable" << std::endl;
            pk = prev.pk + alpha_backstep*vp;
            setParams(pk);
            traj = Trajectory(prev_zk_state);
            alpha_backstep/=2.0;    
        }
        traj = Trajectory({zk.segment(0, NUM_PLANT_STATES), zk.segment(NUM_PLANT_STATES, NUM_ALGE_STATES)});
        backtrackingCount ++;
    }
    if(backtrackingCount > 0)
    {
        std::chrono::time_point end = std::chrono::steady_clock::now();
        std::chrono::microseconds elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
        m_logger << "backtracking time: " << elapsed.count() << " us" << std::endl;
        m_logger << "backtracking count: " << backtrackingCount << std::endl;
        
    }
    if(!traj.stable || !traj.converged){
        return prev;
    } else {
        m_logger << "stable" << std::endl;
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
        while(!traj.stable && (zk-zk_prev).cwiseAbs().sum() > 1e-8)
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