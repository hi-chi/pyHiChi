#include "PicParticleReader.h"
#include <fstream>

void PicParticleReader::ReadFromFile(const std::vector<std::string>& fileNames,
    PicParticleContainer& particles)
{
    size_t counter = 1;
    std::for_each(fileNames.begin(), fileNames.end(),
        [&counter, &particles](const std::string& fileName) {
            // считываем из файла частицы в массив структур PicParticle
            std::ifstream file(fileName);

            if (file.is_open()) {

                if (!file.eof()) { // заголовок
                    std::string line;
                    std::getline(file, line);
                }

                while (!file.eof()) {

                    std::stringstream sstr;
                    std::string line;
                    std::getline(file, line);
                    sstr << line;

                    if (line == "\0" || line == "\n") break;

                    PicParticle particle;
                    particle.id = counter++;

                    std::string stype;
                    sstr >> stype;
                    particle.SetType(stype);

                    sstr >> particle.mass;
                    sstr >> particle.charge;
                    sstr >> particle.factor;
                    
                    sstr >> particle.position.x;
                    sstr >> particle.position.y;
                    sstr >> particle.position.z;
                    
                    sstr >> particle.velocity.x;
                    sstr >> particle.velocity.y;
                    sstr >> particle.velocity.z;
                     
                    sstr >> particle.momentum.x;
                    sstr >> particle.momentum.y;
                    sstr >> particle.momentum.z;
                    
                    sstr >> particle.emissionTime;

                    particles.Push(particle);
                }

                file.close();
            }
            else {
                std::cout << "ERROR: can not open file " << fileName << "." << std::endl;
            }
        });
}