#pragma once
#include "Constants.h"
#include "Ensemble.h"
#include "Grid.h"
#include "AnalyticalField.h"
#include "Pusher.h"
#include "Vectors.h"

#include "compton.h"
#include "breit_wheeler.h"

#include <stdint.h>
#include <omp.h>
#include <random>

using namespace constants;
namespace pfc
{
    template <class TGrid>  // may be AnalyticalField or any Grid type
    class Scalar_Fast_QED : public ParticlePusher
    {
    public:

        Scalar_Fast_QED() : compton(), breit_wheeler()
        {
            SchwingerField = sqr(Constants<FP>::electronMass() * Constants<FP>::lightVelocity())
                * Constants<FP>::lightVelocity() / (-Constants<FP>::electronCharge() * Constants<FP>::planck());

            coeffPhoton_probability = 1.0;
            coeffPair_probability = 1.0;

            distribution = std::uniform_real_distribution<FP>(0.0, 1.0);
            int max_threads;

#ifdef __USE_OMP__
            max_threads = omp_get_max_threads();
#else
            max_threads = 1;
#endif
            AvalanchePhotons.resize(max_threads);
            AvalancheParticles.resize(max_threads);
            afterAvalanchePhotons.resize(max_threads);
            afterAvalancheParticles.resize(max_threads);
            timeAvalanchePhotons.resize(max_threads);
            timeAvalancheParticles.resize(max_threads);
        }

        void disable_photon_emission() 
        {
            this->coeffPhoton_probability = 0.0;
        }

        void enable_photon_emission()
        {
            this->coeffPhoton_probability = 1.0;
        }

        void disable_pair_production()
        {
            this->coeffPair_probability = 0.0;
        }

        void enable_pair_production()
        {
            this->coeffPair_probability = 1.0;
        }

        void processParticles(Ensemble3d* particles, TGrid* grid, FP timeStep)
        {
            int max_threads;
#ifdef __USE_OMP__
            max_threads = omp_get_max_threads();
#else
            max_threads = 1;
#endif

            for (int th = 0; th < max_threads; th++)
            {
                AvalanchePhotons[th].clear();
                AvalancheParticles[th].clear();
                afterAvalanchePhotons[th].clear();
                afterAvalancheParticles[th].clear();
            }

            if ((*particles)[Photon].size())
                HandlePhotons((*particles)[Photon], grid, timeStep);
            if ((*particles)[Electron].size())
                HandleParticles((*particles)[Electron], grid, timeStep);
            if ((*particles)[Positron].size())
                HandleParticles((*particles)[Positron], grid, timeStep);

            for (int th = 0; th < max_threads; th++)
            {
                for (int ind = 0; ind < afterAvalanchePhotons[th].size(); ind++)
                {
                    particles->addParticle(afterAvalanchePhotons[th][ind]);
                }
                for (int ind = 0; ind < afterAvalancheParticles[th].size(); ind++)
                {
                    particles->addParticle(afterAvalancheParticles[th][ind]);
                }
            }
        }

        void Boris(Particle3d&& particle, const FP3& e, const FP3& b, FP timeStep)
        {
            FP eCoeff = timeStep * particle.getCharge() / (2 * particle.getMass() * Constants<FP>::lightVelocity());
            FP3 eMomentum = e * eCoeff;
            FP3 um = particle.getP() + eMomentum;
            FP3 t = b * eCoeff / sqrt((FP)1 + um.norm2());
            FP3 uprime = um + cross(um, t);
            FP3 s = t * (FP)2 / ((FP)1 + t.norm2());
            particle.setP(eMomentum + um + cross(uprime, s));
            particle.setPosition(particle.getPosition() + timeStep * particle.getVelocity());
        }

        void Boris(ParticleProxy3d&& particle, const FP3& e, const FP3& b, FP timeStep)
        {
            FP eCoeff = timeStep * particle.getCharge() / (2 * particle.getMass() * Constants<FP>::lightVelocity());
            FP3 eMomentum = e * eCoeff;
            FP3 um = particle.getP() + eMomentum;
            FP3 t = b * eCoeff / sqrt((FP)1 + um.norm2());
            FP3 uprime = um + cross(um, t);
            FP3 s = t * (FP)2 / ((FP)1 + t.norm2());
            particle.setP(eMomentum + um + cross(uprime, s));
            particle.setPosition(particle.getPosition() + timeStep * particle.getVelocity());
        }

        void HandlePhotons(ParticleArray3d& particles, TGrid* grid, FP timeStep)
        {
#pragma omp parallel for
                for (int i = 0; i < particles.size(); i++)
                {
                    int thread_id;
#ifdef __USE_OMP__
                    thread_id = omp_get_thread_num();
#else
                    thread_id = 0;
#endif
                    FP3 pPos = particles[i].getPosition();
                    FP3 e, b;

                    e = grid->getE(pPos);
                    b = grid->getB(pPos);

                    AvalancheParticles[thread_id].clear();
                    AvalanchePhotons[thread_id].clear();
                    timeAvalancheParticles[thread_id].clear();
                    timeAvalanchePhotons[thread_id].clear();
                    AvalanchePhotons[thread_id].push_back(particles[i]);
                    timeAvalanchePhotons[thread_id].push_back((FP)0.0);

                    RunAvalanche(e, b, timeStep);

                    for (int k = 0; k != AvalanchePhotons[thread_id].size(); k++)
                        if (AvalanchePhotons[thread_id][k].getGamma() != (FP)1.0)
                            afterAvalanchePhotons[thread_id].push_back(AvalanchePhotons[thread_id][k]);


                    for (int k = 0; k != AvalancheParticles[thread_id].size(); k++)
                        afterAvalancheParticles[thread_id].push_back(AvalancheParticles[thread_id][k]);
                }

                particles.clear();
        }

        void HandleParticles(ParticleArray3d& particles, TGrid* grid, FP timeStep)
        {
            
#pragma omp parallel for
            for (int i = 0; i < particles.size(); i++)
            {
                int thread_id;
#ifdef __USE_OMP__
                thread_id = omp_get_thread_num();
#else
                thread_id = 0;
#endif
                FP3 pPos = particles[i].getPosition();
                FP3 e, b;

                e = grid->getE(pPos);
                b = grid->getB(pPos);
                
                AvalancheParticles[thread_id].clear();
                AvalanchePhotons[thread_id].clear();
                timeAvalancheParticles[thread_id].clear();
                timeAvalanchePhotons[thread_id].clear();
                AvalancheParticles[thread_id].push_back(particles[i]);
                timeAvalancheParticles[thread_id].push_back((FP)0.0);

                RunAvalanche(e, b, timeStep);

                for (int k = 0; k != AvalanchePhotons[thread_id].size(); k++)
                    if(AvalanchePhotons[thread_id][k].getGamma() != (FP)1.0)
                        afterAvalanchePhotons[thread_id].push_back(AvalanchePhotons[thread_id][k]);

                particles[i].setMomentum(AvalancheParticles[thread_id][0].getMomentum());
                particles[i].setPosition(AvalancheParticles[thread_id][0].getPosition());

                for (int k = 1; k != AvalancheParticles[thread_id].size(); k++)
                    afterAvalancheParticles[thread_id].push_back(AvalancheParticles[thread_id][k]);
            }
        }
        
        void RunAvalanche(const FP3& E, const FP3& B, double timeStep)
        {
            int thread_id;
#ifdef __USE_OMP__
            thread_id = omp_get_thread_num();
#else
            thread_id = 0;
#endif
            vector<Particle3d>& AvalancheParticles = this->AvalancheParticles[thread_id];
            vector<Particle3d>& AvalanchePhotons = this->AvalanchePhotons[thread_id];

            vector<FP>& timeAvalancheParticles = this->timeAvalancheParticles[thread_id];
            vector<FP>& timeAvalanchePhotons = this->timeAvalanchePhotons[thread_id];
            
            int countParticles = 0;
            int countPhotons = 0;

            while (countParticles != AvalancheParticles.size()
                || countPhotons != AvalanchePhotons.size())
            {
                for (int k = countParticles; k != AvalancheParticles.size(); k++)
                {
                    oneParticleStep(AvalancheParticles[k], E, B, timeAvalancheParticles[k], timeStep);
                    countParticles++;
                }
                for (int k = countPhotons; k != AvalanchePhotons.size(); k++)
                {
                    onePhotonStep(AvalanchePhotons[k], E, B, timeAvalanchePhotons[k], timeStep);
                    countPhotons++;
                }
            }
        }

        void oneParticleStep(Particle3d& particle, const FP3& E, const FP3& B, FP& time, double timeStep)
        {
            int thread_id;
#ifdef __USE_OMP__
            thread_id = omp_get_thread_num();
#else
            thread_id = 0;
#endif
            while (time < timeStep)
            {
                FP3 v = particle.getVelocity();
                FP H_eff = sqr(E + (1 / Constants<FP>::lightVelocity()) * VP(v, B))
                    - sqr(SP(E, v) / Constants<FP>::lightVelocity());
                if (H_eff < 0) H_eff = 0;
                H_eff = sqrt(H_eff);
                FP gamma = particle.getGamma();
                FP chi = gamma * H_eff / SchwingerField;
                FP rate = 0.0, dt = 2*timeStep;
                if (chi > 0.0 && coeffPhoton_probability != 0.0)
                {
                    rate = compton.rate(chi);
                    dt = getDtParticle(particle, rate, chi);
                }

                if (dt + time > timeStep)
                {
                    Boris(particle, E, B, timeStep - time);
                    time = timeStep;
                }
                else
                {
                    Boris(particle, E, B, dt);
                    time += dt;
                    FP3 v = particle.getVelocity();
                    FP H_eff = sqr(E + (1 / Constants<FP>::lightVelocity()) * VP(v, B))
                        - sqr(SP(E, v) / Constants<FP>::lightVelocity());
                    if (H_eff < 0) H_eff = 0;
                    H_eff = sqrt(H_eff);
                    FP gamma = particle.getGamma();
                    FP chi_new = gamma * H_eff / SchwingerField;

                    FP delta = Photon_Generator((chi + chi_new)/(FP)2.0);
                    
                    Particle3d NewParticle;
                    NewParticle.setType(Photon);
                    NewParticle.setWeight(particle.getWeight());
                    NewParticle.setPosition(particle.getPosition());
                    NewParticle.setMomentum(delta * particle.getMomentum());

                    this->AvalanchePhotons[thread_id].push_back(NewParticle);
                    this->timeAvalanchePhotons[thread_id].push_back(time);
                    particle.setMomentum((1 - delta) * particle.getMomentum());
                }
            }
        }

        void onePhotonStep(Particle3d& particle, const FP3& E, const FP3& B, FP& time, double timeStep)
        {
            int thread_id;
#ifdef __USE_OMP__
            thread_id = omp_get_thread_num();
#else
            thread_id = 0;
#endif

            FP3 k_ = particle.getVelocity();
            k_ = (1 / k_.norm()) * k_; // normalized wave vector
            FP H_eff = sqrt(sqr(E + VP(k_, B)) - sqr(SP(E, k_)));
            FP gamma = particle.getMomentum().norm()
                / (Constants<FP>::electronMass() * Constants<FP>::lightVelocity());
            FP chi = gamma * H_eff / SchwingerField;

            FP rate = 0.0, dt = 2 * timeStep;
            if (chi > 0.0 && coeffPair_probability != 0.0)
            {
                rate = breit_wheeler.rate(chi);
                dt = getDtPhoton(particle, rate, chi, gamma);

            }

            if (dt + time > timeStep)
            {
                particle.setPosition(particle.getPosition()
                    + (timeStep - time) * Constants<FP>::lightVelocity() * k_);
                time = timeStep;
            }
            else
            {
                particle.setPosition(particle.getPosition()
                    + dt * Constants<FP>::lightVelocity() * k_);
                time += dt;
                FP delta = Pair_Generator(chi);

                Particle3d NewParticle;
                NewParticle.setType(Electron);
                NewParticle.setWeight(particle.getWeight());
                NewParticle.setPosition(particle.getPosition());
                NewParticle.setMomentum(delta * particle.getMomentum());

                this->AvalancheParticles[thread_id].push_back(NewParticle);
                this->timeAvalancheParticles[thread_id].push_back(time);

                NewParticle.setType(Positron);
                NewParticle.setMomentum((1 - delta) * particle.getMomentum());
                this->AvalancheParticles[thread_id].push_back(NewParticle);
                this->timeAvalancheParticles[thread_id].push_back(time);

                particle.setP((FP)0.0 * particle.getP());
                time = timeStep;
            }
        }
        
        FP getDtParticle(Particle3d& particle, FP rate, FP chi)
        {
            FP r = -log(random_number_omp());
            r *= (particle.getGamma()) / (chi);
            return r / rate;
        }
        
        FP Photon_Generator(FP chi)
        {
            FP r = random_number_omp();
            return compton.inv_cdf(r, chi);
        }


        FP getDtPhoton(Particle3d& particle, FP rate, FP chi, FP gamma)
        {
            FP r = -log(random_number_omp());
            r *= (gamma) / (chi);
            return r / rate;
        }
      
        FP Pair_Generator(FP chi)
        {
            FP r = random_number_omp();
            if (r < 0.5)
                return breit_wheeler.inv_cdf(r, chi);
            else
                return 1.0 - breit_wheeler.inv_cdf(1.0 - r, chi);
        }

        void operator()(ParticleProxy3d* particle, ValueField field, FP timeStep)
        {}

        void operator()(Particle3d* particle, ValueField field, FP timeStep)
        {
            ParticleProxy3d particleProxy(*particle);
            this->operator()(&particleProxy, field, timeStep);
        }


        Compton compton;
        Breit_wheeler breit_wheeler;

    private:

        FP random_number_omp()
        {
            FP rand_n;
#pragma omp critical
            rand_n = distribution(rand_generator);
            return rand_n;
        }

        std::default_random_engine rand_generator;
        std::uniform_real_distribution<FP> distribution;

        vector<vector<FP>> timeAvalanchePhotons, timeAvalancheParticles;
        vector<vector<Particle3d>> AvalanchePhotons, AvalancheParticles;
        vector<vector<Particle3d>> afterAvalanchePhotons, afterAvalancheParticles;

        FP SchwingerField;
        
        FP coeffPhoton_probability, coeffPair_probability;
    };

    typedef Scalar_Fast_QED<YeeGrid> Scalar_Fast_QED_Yee;
    typedef Scalar_Fast_QED<PSTDGrid> Scalar_Fast_QED_PSTD;
    typedef Scalar_Fast_QED<PSATDGrid> Scalar_Fast_QED_PSATD;
    typedef Scalar_Fast_QED<AnalyticalField> Scalar_Fast_QED_Analytical;
}
