#pragma once
#include "Constants.h"
#include "FP.h"
#include "Particle.h"
#include "ParticleArray.h"

#include <algorithm>
#include <random>
#include <vector>

namespace pfc
{
    template<class ParticleArray>
    class Merging
    {
    public:
        static void merge_with_kmeans(ParticleArray& particles, int numClusters, int iteration)
        {
            std::random_device rd;
            std::mt19937 generator(rd());

            //Copy particle momentum
            vector<FP3> momentums(particles.size());
            for (int idx = 0; idx < particles.size(); idx++)
                momentums[idx] = particles[idx].getMomentum();

            // k-means for getting cluster numbers
            vector<int> clusterDecomposition = kMeans(momentums, numClusters, iteration);

            vector<vector<Particle3d>> clusters(numClusters);
            for (int i = 0; i < particles.size(); i++)
                clusters[clusterDecomposition[i]].push_back(Particle3d(particles[i]));

            // Instead of each cluster, we create one particle with a total factor;
            // the position ones from all, momentum are weighted averages, taking into account the factors
            particles.clear();
            for (int j = 0; j < numClusters; j++) {
                if (clusters[j].empty())
                    continue;
                std::uniform_int_distribution<unsigned> uniform_dist(0, clusters[j].size() - 1);
                FP3 position, momentum;
                vector<FP3> positions;
                double weight = 0.0;
                for (int k = 0; k < clusters[j].size(); k++) {
                    positions.push_back(clusters[j][k].getPosition());
                    momentum += clusters[j][k].getMomentum() * clusters[j][k].getWeight();
                    weight += clusters[j][k].getWeight();
                }
                Particle3d mergedParticle(positions[uniform_dist(generator)], momentum / weight,
                    weight, clusters[j][0].getType());
                particles.pushBack(mergedParticle);
            }
        }

    private:
        static vector<int> kMeans(const vector<FP3>& Particles, int k, int iteration)
        {
            std::random_device rd;
            std::mt19937 generator(rd());
            std::uniform_real_distribution<double> uniform_dist(0.0, 1.0);
            int n = (int)Particles.size();
            vector<int> Clusters(n);
            /**
            n is a number of particles
            k is a number of clusters
            Clusters is a 1-dimensional massive containing number of a cluster which contains particle number i
            initial centroids are choosen as random particles from Particles[]
            **/
            int numCyc = iteration; ///maximal number of iterations to determine centroids
            vector<FP3> Center(k);///current coordinates of centroids
            for (int i = 0; i < n; i++)
                Clusters[i] = 1;///temporary initialization for choosing of initial centroids
            ///selection of centroids' coordinates
            Center[0] = Particles[0];
            Clusters[0] = 0;///0 for massive Clusters  means that particle currentNum is used as center of a cluster
            vector<double> minDist(n);///temporary massive containing distance from particle to the nearest centroid
            for (int j = 1; j < k; j++)
            {
                double sum = 0;/// sum of distances from particles to the nearest centroids
                for (int i = 1; i < n; i++)
                {
                    if (Clusters[i] != 0)
                    {
                        minDist[i] = dist(Particles[i], Center[0]);
                        double buf = 0;
                        for (int h = 1; h < j; h++)
                        {
                            buf = dist(Particles[i], Center[h]);
                            if (buf < minDist[i])
                            {
                                minDist[i] = buf;
                            }
                        }
                        sum += minDist[i];
                    }
                }

                double threshChoise = uniform_dist(generator) * sum;
                double tsum = 0;
                int np = 1;
                while (tsum < threshChoise)
                {
                    if (Clusters[np] != 0)
                    {
                        tsum += minDist[np];
                    }
                    np++;
                }
                Center[j] = Particles[np - 1];
                Clusters[np - 1] = 0;
            }

            ///main loop, if this algorithm doesn't converge in numCyc iterations then it is interrupted
            for (int j = 0; j < numCyc; j++)
            {
                vector<int> nPart(k, 0);///massive contains quantity of particles belonging to the cluster
                vector<double> sumDist(k, 0.0);///for estimation of dispersion
                vector<FP3> newCenter(k);///creation variables for calculating new coordinates of centroids
                for (int np = 0; np < n; np++) ///loop for determining particles of clusters
                {
                    int numcl = 0;///number of a cluster contains particle number np
                    double d = dist(Particles[np], Center[0]);
                    for (int nc = 1; nc < k; nc++) ///search for nearest centroid
                    {
                        double buf = dist(Particles[np], Center[nc]);
                        if (d > buf)
                        {
                            d = buf;
                            numcl = nc;
                        }
                    }
                    Clusters[np] = numcl;
                    nPart[numcl] += 1;
                    newCenter[numcl] += Particles[np];
                    sumDist[numcl] += d;
                }
                bool isSame = true;
                double max = 0;
                for (int nc = 0; nc < k; nc++)
                {
                    if (nPart[nc] > 0)
                    {
                        newCenter[nc] = newCenter[nc] / static_cast<double>(nPart[nc]);
                        isSame &= (newCenter[nc] == Center[nc]);
                        Center[nc] = newCenter[nc];
                        sumDist[nc] = sqrt(sumDist[nc] / static_cast<double>(nPart[nc]));
                        if (max < sumDist[nc])
                            max = sumDist[nc];
                    }

                }
                if (isSame)
                    break;
            }
            return Clusters;
        }
    };
}