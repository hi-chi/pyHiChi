#pragma once
#include "Constants.h"
#include "FP.h"
#include "ParticleArray.h"

#include <algorithm>
#include <random>
#include <vector>
#include <queue>

namespace pfc
{
    template<class ParticleArray>
    class Thinning
    {

        std::random_device rd;
        std::mt19937 generator;
    public:
        Thinning() : generator(rd())
        {}

        void simple(ParticleArray& particles, int m)
        {
            int sizeArray = particles.size();

            for (int i = 0; i < sizeArray - m; i++)
            {
                std::uniform_int_distribution<int> dist(0, sizeArray - i - 1);
                int deleteIdx = dist(generator);
                particles.deleteParticle(deleteIdx);
            }
            FP newCoeff = static_cast<FP>(sizeArray) / (static_cast<FP>(m));
            for (int idx = 0; idx < particles.size(); idx++)
            {
                FP weight = particles[idx].getWeight();
                particles[idx].setWeight(weight * newCoeff);
            }
        }

        void leveling(ParticleArray& particles)
        {
            std::uniform_real_distribution<FP> dist(0, 1);

            FP weightAvg = 0.0;
            for (int idx = 0; idx < particles.size(); idx++)
            {
                weightAvg += particles[idx].getWeight() / particles.size();
            }
            FP threshold = 2 * weightAvg;

            for (int idx = particles.size() - 1; idx >= 0; idx--)
            {
                if (particles[idx].getWeight() < threshold)
                {
                    FP randNumber = dist(generator);
                    if (randNumber > particles[idx].getWeight() / threshold)
                    {
                        particles.deleteParticle(idx);
                    }
                    else
                    {
                        particles[idx].setWeight(threshold);
                    }
                }
            }
        }

        void numberConservative(ParticleArray& particles, int m)
        {
            FP weightSum = 0.0;
            for (int idx = 0; idx < particles.size(); idx++)
            {
                weightSum += particles[idx].getWeight();
            }
            std::uniform_real_distribution<FP> dist(0, weightSum);
            std::vector<FP> randomNumbers;
            randomNumbers.resize(m + 1);
            for (int idx = 0; idx < m; idx++)
            {
                randomNumbers[idx] = dist(generator);
            }
            randomNumbers[m] = weightSum * 2;
            std::sort(randomNumbers.begin(), randomNumbers.end());

            int idxNumber = 0;
            FP currentWeightSum = 0.0;
            for (int idx = particles.size() - 1; idx >= 0; idx--)
            {
                int ki = 0;
                currentWeightSum += particles[idx].getWeight();
                while (currentWeightSum > randomNumbers[idxNumber])
                {
                    ki++;
                    idxNumber++;
                }
                if (ki == 0)
                {
                    particles.deleteParticle(idx);
                }
                else
                {
                    FP newCoeff = static_cast<FP>(ki) / static_cast<FP>(m);
                    particles[idx].setWeight(newCoeff * weightSum);
                }
            }
        }

        void energyConservative(ParticleArray& particles, int m)
        {
            FP energySum = 0.0;
            std::vector<FP> energys;
            for (int idx = 0; idx < particles.size(); idx++)
            {
                FP mass = particles[idx].getMass();
                FP c = Constants<FP>::lightVelocity();
                FP mc = mass * c;
                FP energyParticle = sqrt(mc * mc + particles[idx].getMomentum().norm2());
                energys.push_back(energyParticle);
                energySum += particles[idx].getWeight() * energyParticle;
            }
            std::uniform_real_distribution<FP> dist(0, energySum);
            std::vector<FP> randomNumbers;
            randomNumbers.resize(m + 1);
            for (int idx = 0; idx < m; idx++)
            {
                randomNumbers[idx] = dist(generator);
            }
            randomNumbers[m] = energySum * 2;
            std::sort(randomNumbers.begin(), randomNumbers.end());

            int idxNumber = 0;
            FP currentEnergySum = 0.0;
            for (int idx = particles.size() - 1; idx >= 0; idx--)
            {
                int ki = 0;
                currentEnergySum += particles[idx].getWeight() * energys[idx];
                while (currentEnergySum > randomNumbers[idxNumber])
                {
                    ki++;
                    idxNumber++;
                }
                if (ki == 0)
                {
                    particles.deleteParticle(idx);
                }
                else
                {
                    FP newCoeff = static_cast<FP>(ki) / static_cast<FP>(m);
                    particles[idx].setWeight(newCoeff * energySum / energys[idx]);
                }
            }
        }

        enum Features
        {
            Energy,
            Momentum, Momentum_x, Momentum_y, Momentum_z,
            Position, Position_x, Position_y, Position_z,
            Dispersion_Energy,
            Dispersion_Momentum, Dispersion_Momentum_x, 
            Dispersion_Momentum_y, Dispersion_Momentum_z,
            Dispersion_Position, Dispersion_Position_x,
            Dispersion_Position_y, Dispersion_Position_z,
        };

        void thinningConservative(ParticleArray& particles, int m, std::set<Features>& features)
        {
            auto searchM = features.find(Momentum);
            auto searchP = features.find(Position);
            auto searchDM = features.find(Dispersion_Momentum);
            auto searchDP = features.find(Dispersion_Position);
            if (searchM != features.end()) {
                features.erase(Momentum);
                features.insert(Momentum_x);
                features.insert(Momentum_y);
                features.insert(Momentum_z);
            }
            if (searchP != features.end()) {
                features.erase(Position);
                features.insert(Position_x);
                features.insert(Position_y);
                features.insert(Position_z);
            }
            if (searchDM != features.end()) {
                features.erase(Dispersion_Momentum);
                features.insert(Dispersion_Momentum_x);
                features.insert(Dispersion_Momentum_y);
                features.insert(Dispersion_Momentum_z);
            }
            if (searchDP != features.end()) {
                features.erase(Dispersion_Position);
                features.insert(Dispersion_Position_x);
                features.insert(Dispersion_Position_y);
                features.insert(Dispersion_Position_z);
            }

            std::vector<double> weight;
            std::vector< vector<double> > components;
            std::vector<Particle3d> partParticleArray;
            weight.resize(features.size() + 2);
            partParticleArray.resize(features.size() + 2);
            components.resize(features.size() + 1);
            for (int indComp = 0; indComp < features.size(); indComp++)
            {
                components[indComp].resize(features.size() + 2);
            }
            vector<double> unit(features.size() + 2, 1);
            components[components.size() - 1] = unit;

            double c = Constants<FP>::lightVelocity();

            struct particleCompare
            {
                bool operator()(const Particle3d& lhs, const Particle3d& rhs)
                {
                    return lhs.getWeight() > rhs.getWeight();
                }
            };
            typedef std::priority_queue<Particle3d, vector<Particle3d>, particleCompare> queueParticle;
            queueParticle particles_arrays;


            FP mean_energy = (FP)0.0;
            FP3 mean_position = FP3(0.0, 0.0, 0.0);
            FP3  mean_momentum = FP3(0.0, 0.0, 0.0);

            for (int particleIdx = 0; particleIdx < particles.size(); particleIdx++) {

                particles_arrays.push(particles[particleIdx]);

                mean_energy += particles[particleIdx].getWeight() * sqrt(c * c *
                    particles[particleIdx].getMass() * particles[particleIdx].getMass()
                    + particles[particleIdx].getMomentum().norm2());
                mean_position += particles[particleIdx].getWeight() * particles[particleIdx].getPosition();
                mean_momentum += particles[particleIdx].getWeight() * particles[particleIdx].getMomentum();

                particles.deleteParticle(particleIdx);
                particleIdx--;
            }
            mean_energy /= (FP)particles_arrays.size();
            mean_position /= (FP)particles_arrays.size();
            mean_momentum /= (FP)particles_arrays.size();

            int reqDel = particles_arrays.size() * (1.0 - (FP)m / (FP)particles_arrays.size());
            int delParticles = 0;
            while ((particles_arrays.size() >= partParticleArray.size()) &&
                (delParticles < reqDel))
            {
                double mass = particles_arrays.top().getMass();
                double mc2 = mass * c * mass * c;
                for (int indexPart = 0; indexPart < partParticleArray.size(); indexPart++)
                {
                    partParticleArray[indexPart] = particles_arrays.top();
                    particles_arrays.pop();
                    weight[indexPart] = partParticleArray[indexPart].getWeight();
                    int indexFeat = 0;
                    for (auto it = features.begin(); it != features.end();
                        ++it, indexFeat++)
                    {
                        switch (*it)
                        {
                        case Energy:
                            components[indexFeat][indexPart] = sqrt(mc2 + partParticleArray[indexPart].getMomentum().norm2());
                            break;
                        case Momentum_x:
                            components[indexFeat][indexPart] = partParticleArray[indexPart].getMomentum().x;
                            break;
                        case Momentum_y:
                            components[indexFeat][indexPart] = partParticleArray[indexPart].getMomentum().y;
                            break;
                        case Momentum_z:
                            components[indexFeat][indexPart] = partParticleArray[indexPart].getMomentum().z;
                            break;
                        case Position_x:
                            components[indexFeat][indexPart] = partParticleArray[indexPart].getPosition().x;
                            break;
                        case Position_y:
                            components[indexFeat][indexPart] = partParticleArray[indexPart].getPosition().y;
                            break;
                        case Position_z:
                            components[indexFeat][indexPart] = partParticleArray[indexPart].getPosition().z;
                            break;

                        case Dispersion_Energy:
                            components[indexFeat][indexPart] =
                                sqr(sqrt(mc2 + partParticleArray[indexPart].getMomentum().norm2()) - mean_energy);
                            break;
                        case Dispersion_Momentum_x:
                            components[indexFeat][indexPart] =
                                sqr(partParticleArray[indexPart].getMomentum().x - mean_momentum.x);
                            break;
                        case Dispersion_Momentum_y:
                            components[indexFeat][indexPart] =
                                sqr(partParticleArray[indexPart].getMomentum().y - mean_momentum.y);
                            break;
                        case Dispersion_Momentum_z:
                            components[indexFeat][indexPart] =
                                sqr(partParticleArray[indexPart].getMomentum().z - mean_momentum.z);
                            break;
                        case Dispersion_Position_x:
                            components[indexFeat][indexPart] =
                                sqr(partParticleArray[indexPart].getPosition().x - mean_position.x);
                            break;
                        case Dispersion_Position_y:
                            components[indexFeat][indexPart] =
                                sqr(partParticleArray[indexPart].getPosition().y - mean_position.y);
                            break;
                        case Dispersion_Position_z:
                            components[indexFeat][indexPart] =
                                sqr(partParticleArray[indexPart].getPosition().z - mean_position.z);
                            break;
                        }
                    }
                }
                //checked system
                std::vector<double> weight1;
                std::vector< vector<double> > components1;
                for (int indComp = 0; indComp < components.size() - 1; indComp++)
                {
                    double max_elem = std::abs(components[indComp][0]);
                    double min_elem = std::abs(components[indComp][0]);
                    for (int indexPart = 1; indexPart < partParticleArray.size(); indexPart++)
                    {
                        if (max_elem < std::abs(components[indComp][indexPart]))
                            max_elem = std::abs(components[indComp][indexPart]);
                        else if (min_elem > std::abs(components[indComp][indexPart]))
                            min_elem = std::abs(components[indComp][indexPart]);

                    }
                    if (min_elem != max_elem)
                    {
                        for (int indexPart = 0; indexPart < partParticleArray.size(); indexPart++)
                            components[indComp][indexPart] /= max_elem;
                        components1.push_back(components[indComp]);
                    }

                }

                components1.push_back(components[components.size() - 1]); //unit

                for (int indPart = 0; indPart <= components1.size(); indPart++)
                {
                    weight1.push_back(weight[indPart]);
                }


                int isNoThin = thinningConservative_reweighing(weight1, components1);
                if (isNoThin == 0)
                {
                    for (int indexPart = 0; indexPart < partParticleArray.size(); indexPart++)
                    {
                        if (indexPart < weight1.size())
                        {
                            if (weight1[indexPart] != 0.0)
                            {
                                partParticleArray[indexPart].setWeight(weight1[indexPart]);
                                particles_arrays.push(partParticleArray[indexPart]);
                            }
                        }
                        else
                        {
                            partParticleArray[indexPart].setWeight(weight[indexPart]);
                            particles_arrays.push(partParticleArray[indexPart]);
                        }
                    }
                    delParticles++;
                }
                else
                {
                    int indexPart;
                    for (indexPart = 0; indexPart < partParticleArray.size() - 1; indexPart++)
                    {
                        particles_arrays.push(partParticleArray[indexPart]);
                    }
                    particles.pushBack(partParticleArray[indexPart]);
                }
            }

            while (particles_arrays.size())
            {
                particles.pushBack(particles_arrays.top());
                particles_arrays.pop();
            }
        }

        inline double dot(vector<double> v1, vector<double> v2)
        {
            double S = 0;
            for (int i = 0; i < v1.size(); i++)
                S += v1[i] * v2[i];
            return S;
        }
        inline vector<double> mult(double multiplier, vector<double> v)
        {
            vector<double> R(v);
            for (int i = 0; i < R.size(); i++)R[i] *= multiplier;
            return R;
        }
        inline vector<double> sum(vector<double> v1, vector<double> v2)
        {
            vector<double> R(v1);
            for (int i = 0; i < R.size(); i++)R[i] = v1[i] + v2[i];
            return R;
        }

        inline int thinningConservative_reweighing(vector<double>& weight, vector< vector<double> >& components)
        {
            std::uniform_real_distribution<FP> dist(0, 1.0);
         
            int N = weight.size();

            vector< vector<double> > v(N);
            for (int i = 0; i < N; i++)
            {
                v[i].assign(N, 0);
                v[i][i] = 1;
            }

            for (int inv = 0; inv < components.size(); inv++)
            {
                for (int j = 0; j < v.size() - 1; j++)
                {
                    double k1 = dot(v[j + 1], components[inv]);
                    double k2 = -dot(v[j], components[inv]);
                    if ((k1 != 0.0) || (k2 != 0.0))
                        v[j] = sum(mult(k1, v[j]), mult(k2, v[j + 1]));
                    else
                        return -1;

                    double norm = 0;
                    for (int ind = 0; ind < v[j].size(); ind++)
                        norm += v[j][ind] * v[j][ind];
                    norm = sqrt(norm);
                    if (norm != 0.0)
                        for (int ind = 0; ind < v[j].size(); ind++)
                            v[j][ind] /= norm;
                    else
                        return -2;
                }
                v.pop_back();
            }

            // simple choice
            vector<double> r(v[0]);

            double ap, an; //positive and negrative values of a that are closest to 0;
            int ip, in; // the corresponding indices;
            bool fp = false, fn = false; // found negative, found positive
            for (int i = 0; i < N; i++)
            {
                if (r[i] != 0.0)
                {
                    double a = -weight[i] / r[i];
                    if (a > 0)
                    {
                        if (!fp)
                        {
                            fp = true;  ap = a;  ip = i;
                        }
                        if (a < ap)
                        {
                            ap = a;  ip = i;
                        }
                    }
                    else
                    {
                        if (!fn)
                        {
                            fn = true;  an = a;  in = i;
                        }
                        if (a > an)
                        {
                            an = a;  in = i;
                        }
                    }
                }
            }
            if (fp && fn)
            {
                if (dist(generator) < ap / (ap - an))
                {
                    weight = sum(weight, mult(an, r));
                    weight[in] = 0;
                }
                else
                {
                    weight = sum(weight, mult(ap, r));
                    weight[ip] = 0;
                }
                return 0;
            }
            else
                return -3;
        }
    };
}
