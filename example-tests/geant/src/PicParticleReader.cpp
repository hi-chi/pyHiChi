#include "PicParticleReader.h"
#include <fstream>

void PicParticleReader::ReadFromFile(const std::vector<std::string>& fileNames,
    PicParticleContainer& particles)
{
    std::for_each(fileNames.begin(), fileNames.end(),
        [&particles](const std::string& fileName) {
            // считываем из файла частицы в массив структур PicParticle
            std::ifstream file(fileName, std::ios::binary);

            if (file.is_open()) {

                std::string stype;
                // read type
                if (!file.eof()) {
                    std::getline(file, stype);
                }
                else throw(std::exception("particle type not found, incorrect format"));

                size_t nParticles = 0;

                try {
                    // read and save particles
                    while (!file.eof()) {

                        PicParticle particle;

                        particle.id = nParticles++;
                        particle.SetType(stype);

                        double mass = 0.0, charge = 0.0, factor = 0.0,
                            posx = 0.0, posy = 0.0, posz = 0.0,
                            velx = 0.0, vely = 0.0, velz = 0.0,
                            momx = 0.0, momy = 0.0, momz = 0.0, t = 0.0;

                        file.read((char*)&mass, sizeof(mass));
                        file.read((char*)&charge, sizeof(charge));
                        file.read((char*)&factor, sizeof(factor));

                        file.read((char*)&posx, sizeof(posx));
                        file.read((char*)&posy, sizeof(posy));
                        file.read((char*)&posz, sizeof(posz));

                        file.read((char*)&velx, sizeof(velx));
                        file.read((char*)&vely, sizeof(vely));
                        file.read((char*)&velz, sizeof(velz));

                        file.read((char*)&momx, sizeof(momx));
                        file.read((char*)&momy, sizeof(momy));
                        file.read((char*)&momz, sizeof(momz));

                        file.read((char*)&t, sizeof(t));

                        particle.mass = mass;
                        particle.charge = charge;
                        particle.factor = factor;

                        particle.position.x = posx;
                        particle.position.y = posy;
                        particle.position.z = posz;

                        particle.velocity.x = velx;
                        particle.velocity.y = vely;
                        particle.velocity.z = velz;

                        particle.momentum.x = momx;
                        particle.momentum.y = momy;
                        particle.momentum.z = momz;

                        particle.emissionTime = t;

                        particles.Push(particle);
                    }
                }
                catch (...) {
                    throw std::exception("error when reading particles");
                }

                file.close();
            }
            else throw std::exception(("can not open file " + fileName).c_str());
        });
}