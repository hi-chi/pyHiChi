#pragma once
#include "Constants.h"
#include "FP.h"
#include "ParticleArray.h"

#include <algorithm>
#include <random>
#include <vector>


template<class ParticleArray>
class Thinning
{
public:
	static void simple(ParticleArray& particles, int m)
	{
		std::random_device rd;
		std::mt19937 generator(rd());
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

	static void leveling(ParticleArray& particles)
	{
		std::random_device rd;
		std::mt19937 generator(rd());
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

	static void numberConservative(ParticleArray& particles, int m)
	{
		std::random_device rd;
		std::mt19937 generator(rd());

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

	static void energyConservative(ParticleArray& particles, int m)
	{
		std::random_device rd;
		std::mt19937 generator(rd());

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
};