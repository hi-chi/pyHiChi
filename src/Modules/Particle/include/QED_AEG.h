#pragma once
#include "Constants.h"
#include "Ensemble.h"
#include "Grid.h"
#include "AnalyticalField.h"
#include "Pusher.h"
#include "synchrotron.h"

#include <omp.h>
#include <random>

using namespace constants;
namespace pfc
{
    template <class TGrid>  // may be AnalyticalField or any Grid type
    class ScalarQED_AEG: public ParticlePusher
    {
    public:

        ScalarQED_AEG()
        {
            MinProbability = 5e-4;
            MaxProbability = 0.01;

            SchwingerField = sqr(Constants<FP>::electronMass() * Constants<FP>::lightVelocity())
                * Constants<FP>::lightVelocity() / (-Constants<FP>::electronCharge() * Constants<FP>::planck());

            preFactor = sqr(Constants<FP>::electronCharge()) * Constants<FP>::electronMass()
                * Constants<FP>::lightVelocity() / sqr(Constants<FP>::planck());

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
            FP dt = timeStep;
#pragma omp parallel for schedule(dynamic, 1)
            for (int i = 0; i < particles.size(); i++)
            {
                int thread_id;
#ifdef __USE_OMP__
                thread_id = omp_get_thread_num();
#else
                thread_id = 0;
#endif
                FP3 pPos = particles[i].getPosition();
                FP3 k = particles[i].getVelocity();
                FP3 e, b;

                e = grid->getE(pPos);
                b = grid->getB(pPos);

                k = (1 / k.norm()) * k; // normalized wave vector
                particles[i].setPosition(pPos + dt * Constants<FP>::lightVelocity() * k);

                FP H_eff = sqrt(sqr(e + VP(k, b)) - sqr(SP(e, k)));

                FP HE = H_eff / SchwingerField;
                FP pGamma = particles[i].getMomentum().norm() / (Constants<FP>::electronMass() * Constants<FP>::lightVelocity());
                FP EstimatedProbability = dt * estimatedPhotons(HE, pGamma);

                FP Factor = 1;
                if (EstimatedProbability < MinProbability)
                {
                    FP r0 = random_number_omp();
                    if (r0 > EstimatedProbability / MinProbability)
                    {
                        afterAvalanchePhotons[thread_id].push_back(particles[i]);
                        continue;
                    }
                    else
                        Factor = MinProbability / EstimatedProbability;
                }
                if (EstimatedProbability < MaxProbability)
                {
                    //=======handle single event========
                    double gamma = pGamma;
                    double chi = gamma * H_eff / SchwingerField;
                    double delta = Pair_Generator(Factor, chi, gamma, dt);
                    if (delta != 0)
                    {
                        Particle3d NewParticle;
                        NewParticle.setType(Electron);
                        NewParticle.setWeight(particles[i].getWeight());
                        NewParticle.setPosition(particles[i].getPosition());
                        NewParticle.setMomentum(delta * particles[i].getMomentum());

                        afterAvalancheParticles[thread_id].push_back(NewParticle);

                        NewParticle.setType(Positron);
                        NewParticle.setMomentum((1 - delta) * particles[i].getMomentum());

                        afterAvalancheParticles[thread_id].push_back(NewParticle);
                    }
                    else
                    {
                        afterAvalanchePhotons[thread_id].push_back(particles[i]);
                    }
                }
                else {
                    //=======handle avalanche========
                    
                    particles[i].setPosition(particles[i].getPosition() - dt * Constants<FP>::lightVelocity() * k); // go back
                    AvalancheParticles[thread_id].clear();
                    AvalanchePhotons[thread_id].clear();
                    AvalanchePhotons[thread_id].push_back(particles[i]);
                    
                    RunAvalanche(H_eff, e, b, Photon, pGamma, dt);
                    
                    for (int k = 0; k != AvalanchePhotons[thread_id].size(); k++)
                        afterAvalanchePhotons[thread_id].push_back(AvalanchePhotons[thread_id][k]);
                    for (int k = 0; k != AvalancheParticles[thread_id].size(); k++)
                        afterAvalancheParticles[thread_id].push_back(AvalancheParticles[thread_id][k]);
                    
                }
            }

            particles.clear();
        }

        void HandleParticles(ParticleArray3d& particles, TGrid* grid, FP timeStep)
        {
            FP dt = timeStep;
#pragma omp parallel for schedule(dynamic, 1)
            for (int i = 0; i < particles.size(); i++)
            {
                int thread_id;
#ifdef __USE_OMP__
                thread_id = omp_get_thread_num();
#else
                thread_id = 0;
#endif
                FP3 pPos = particles[i].getPosition();
                FP3 v = particles[i].getVelocity();
                FP3 e, b;

                e = grid->getE(pPos);
                b = grid->getB(pPos);

                FP H_eff = sqr(e + (1 / Constants<FP>::lightVelocity()) * VP(v, b))
                    - sqr(SP(e, v) / Constants<FP>::lightVelocity());
                if (H_eff < 0)
                    H_eff = 0;
                H_eff = sqrt(H_eff);

                FP pGamma = particles[i].getGamma();
                FP HE = H_eff / SchwingerField;
                FP EstimatedProbability = dt * estimatedParticles(HE, pGamma);

                FP Factor = 1;

                if (EstimatedProbability < MinProbability)
                {
                    FP r0 = random_number_omp();
                    if (r0 > EstimatedProbability / MinProbability)
                    {
                        Boris(particles[i], e, b, dt);
                        continue;
                    }
                    else
                        Factor = MinProbability / EstimatedProbability;
                }
                if (EstimatedProbability < MaxProbability)
                {
                    //=======handle single event========
                    double gamma = pGamma;
                    double chi = gamma * H_eff / SchwingerField;
                    double delta = Photon_MGenerator(Factor, chi, gamma, dt);
                    if (delta != 0)
                    {
                        Particle3d NewParticle;
                        NewParticle.setType(Photon);
                        NewParticle.setWeight(particles[i].getWeight());
                        NewParticle.setPosition(particles[i].getPosition());
                        NewParticle.setMomentum(delta * particles[i].getMomentum());

                        afterAvalanchePhotons[thread_id].push_back(NewParticle);

                        particles[i].setMomentum((1 - delta) * particles[i].getMomentum());
                    }
                    Boris(particles[i], e, b, dt);
                }
                else
                {
                    //=======handle avalanche========
                    AvalancheParticles[thread_id].clear();
                    AvalanchePhotons[thread_id].clear();
                    AvalancheParticles[thread_id].push_back(particles[i]);
                    RunAvalanche(H_eff, e, b, particles[i].getType(), pGamma, dt);

                    for (int k = 0; k != AvalanchePhotons[thread_id].size(); k++)
                        afterAvalanchePhotons[thread_id].push_back(AvalanchePhotons[thread_id][k]);

                    particles[i].setMomentum(AvalancheParticles[thread_id][0].getMomentum());
                    particles[i].setPosition(AvalancheParticles[thread_id][0].getPosition());

                    for (int k = 1; k != AvalancheParticles[thread_id].size(); k++)
                        afterAvalancheParticles[thread_id].push_back(AvalancheParticles[thread_id][k]);
                }
            }
        }

        void RunAvalanche(double H_eff_global, const FP3& E, const FP3& B, int SeedType, double gamma, double dt)
        {
            int thread_id;
#ifdef __USE_OMP__
            thread_id = omp_get_thread_num();
#else
            thread_id = 0;
#endif
            vector<Particle3d>& AvalancheParticles = this->AvalancheParticles[thread_id];
            vector<Particle3d>& AvalanchePhotons = this->AvalanchePhotons[thread_id];
            gamma = max(gamma, 1.0);
            FP HE = H_eff_global / SchwingerField;
            FP sub_dt = MaxProbability / estimatedParticles(HE, gamma);
            int NT = 1 + int(dt / sub_dt);
            sub_dt = dt / FP(NT);

            for (int i = 0; i != NT; i++)
            {
                for (int k = 0; k != AvalancheParticles.size(); k++)
                {
                    Boris(AvalancheParticles[k], E, B, sub_dt);

                    FP3 v = AvalancheParticles[k].getVelocity();
                    FP H_eff = sqr(E + (1 / Constants<FP>::lightVelocity()) * VP(v, B))
                        - sqr(SP(E, v) / Constants<FP>::lightVelocity());
                    if (H_eff < 0) H_eff = 0;
                    H_eff = sqrt(H_eff);
                    FP gamma = AvalancheParticles[k].getGamma();
                    FP chi = gamma * H_eff / SchwingerField;
                    FP delta = Photon_MGenerator(1, chi, gamma, sub_dt);
                    if (delta != 0)
                    {
                        Particle3d NewParticle;
                        NewParticle.setType(Photon);
                        NewParticle.setWeight(AvalancheParticles[k].getWeight());
                        NewParticle.setPosition(AvalancheParticles[k].getPosition());
                        NewParticle.setMomentum(delta * AvalancheParticles[k].getMomentum());

                        AvalanchePhotons.push_back(NewParticle);
                        AvalancheParticles[k].setMomentum((1 - delta) * AvalancheParticles[k].getMomentum());
                    }
                }
                for (int k = 0; k < AvalanchePhotons.size(); k++)
                {
                    FP3 k_ = AvalanchePhotons[k].getVelocity();
                    k_ = (1 / k_.norm()) * k_; // normalized wave vector
                    AvalanchePhotons[k].setPosition(AvalanchePhotons[k].getPosition()
                        + sub_dt * Constants<FP>::lightVelocity() * k_);
                    FP H_eff = sqrt(sqr(E + VP(k_, B)) - sqr(SP(E, k_)));
                    FP gamma = AvalanchePhotons[k].getMomentum().norm()
                        / (Constants<FP>::electronMass() * Constants<FP>::lightVelocity());
                    FP chi = gamma * H_eff / SchwingerField;
                    FP delta = Pair_Generator(1, chi, gamma, sub_dt);
                    if (delta != 0)
                    {
                        Particle3d NewParticle;
                        NewParticle.setType(Electron);
                        NewParticle.setWeight(AvalanchePhotons[k].getWeight());
                        NewParticle.setPosition(AvalanchePhotons[k].getPosition());
                        NewParticle.setMomentum(delta * AvalanchePhotons[k].getMomentum());

                        AvalancheParticles.push_back(NewParticle);

                        NewParticle.setType(Positron);
                        NewParticle.setMomentum((1 - delta) * AvalanchePhotons[k].getMomentum());
                        AvalancheParticles.push_back(NewParticle);
                        AvalanchePhotons[k] = AvalanchePhotons[AvalanchePhotons.size() - 1];
                        AvalanchePhotons.pop_back();
                        k--;
                    }
                }
            }
        }

        FP estimatedPhotons(FP HE, FP gamma)
        {
            return (0.0827 * HE) * preFactor * coeffPair_probability;
        }

        FP estimatedParticles(FP HE, FP gamma)
        {
            FP b = 3.0 / 2.0 * HE * gamma;
            FP newFactor;
            if (b < 0.1)
            {
                newFactor = 0.962436 * b / gamma + 0.0827 * HE;
            }
            else if (b < 0.5)
            {
                newFactor = 0.779009 * pow(b, 11.0 / 12.0) / gamma + 0.0827 * HE;
            }
            else if (b < 10)
            {
                newFactor = 0.721193 * pow(b, 19.0 / 24.0) / gamma + 0.0827 * HE;
            }
            else
            {
                newFactor = 0.955556 * pow(b, 2.0 / 3.0) / gamma + 0.0827 * HE;
            }

            return newFactor * preFactor;
        }


        FP Photon_probability(FP chi, FP gamma, FP d)
        {
            FP z = (2 / 3.0) * (1 / chi) * d / (1 - d);
            FP coeff = (sqrt(3.0) / (2.0 * pi)) * coeffPhoton_probability;
            if ((z < 700) && (z > 0))
                return coeff * (chi / gamma) * ((1 - d) / d) * (synchrotron_1(z) + (3 / 2.0) * d * chi * z * synchrotron_2(z));
            else
                return 0;
        }

        FP Pair_probability(FP chi, FP gamma, FP d)
        {
            FP z_p = (2 / 3.0) / (chi * (1 - d) * d);
            FP coeff = (sqrt(3.0) / (2.0 * pi)) * coeffPair_probability;
            if ((z_p < 700) && (z_p > 0))
                return coeff * (chi / gamma) * (d - 1) * d * (synchrotron_1(z_p) - (3 / 2.0) * chi * z_p * synchrotron_2(z_p));
            else
                return 0;
        }

        FP Pair_Generator(FP Factor, FP chi, FP gamma, FP dt) //returns photon energy in mc2gamma in case of generation.
        {
            FP factor = Factor * dt * preFactor;
            FP r1 = random_number_omp();
            FP r2 = random_number_omp();
            if (r2 < factor * Pair_probability(chi, gamma, r1))
                return r1;
            else
                return 0;
        }
        FP Photon_MGenerator(FP Factor, FP chi, FP gamma, FP dt) //Modified event generator: returns photon energy in mc2gamma in case of generation, !doesn't change gamma
        {
            double r0 = random_number_omp();
            double r1 = r0 * r0 * r0;
            double r2 = random_number_omp();
            double factor = Factor * dt * preFactor;
            if (r2 < factor * Photon_probability(chi, gamma, r1) * 3 * r0 * r0)
                return r1;
            else
                return 0;
        }


        void operator()(ParticleProxy3d* particle, ValueField field, FP timeStep)
        {}

        void operator()(Particle3d* particle, ValueField field, FP timeStep)
        {
            ParticleProxy3d particleProxy(*particle);
            this->operator()(&particleProxy, field, timeStep);
        }
    private:


        FP random_number_omp()
        {
            FP rand_n;
#pragma omp critical
            rand_n = distribution(rand_generator);
            return rand_n;
        }

        FP MinProbability, MaxProbability;
        FP SchwingerField;
        FP preFactor;
        FP coeffPhoton_probability, coeffPair_probability;

        std::default_random_engine rand_generator;
        std::uniform_real_distribution<FP> distribution;


        vector<vector<Particle3d>> AvalanchePhotons, AvalancheParticles;
        vector<vector<Particle3d>> afterAvalanchePhotons, afterAvalancheParticles;
    };

    typedef ScalarQED_AEG<YeeGrid> ScalarQED_AEG_Yee;
    typedef ScalarQED_AEG<PSTDGrid> ScalarQED_AEG_PSTD;
    typedef ScalarQED_AEG<PSATDGrid> ScalarQED_AEG_PSATD;
    typedef ScalarQED_AEG<AnalyticalField> ScalarQED_AEG_Analytical;
}