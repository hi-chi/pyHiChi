#pragma once

#include "ParticleArray.h"
#include <algorithm>
#include <vector>

namespace pfc {
    template <class TGrid, class ArrayType>
    class ParticleSorting {
    public:
        ParticleSorting(TGrid* _grid, ArrayType& particleArray): 
            nn(_grid->numCells.x * _grid->numCells.y * _grid->numCells.z) {
            this->grid = _grid;
            h.resize(nn + 1);
            int idx = 0;
            for (int i = 0; i < particleArray.size(); ++i) {
                FP3 position = particleArray[i].getPosition();
                int linIdx = linIndex(this->grid->getBaseIndex(position));
                while (linIdx >= idx) {
                    h[idx] = i;
                    idx++;
                    
                }
            }
            for (int i = idx; i < h.size(); ++i) {
                h[i] = particleArray.size();
            }
        }

        void sortL(ArrayType& particleArray) {
            for (int n = 0; n < nn - 1; ++n) {
                sortDiap(h, particleArray, n, this->h[n], this->h[n + 1]);
            }   
        }
        
    private:
        TGrid* grid;
        const int nn;
        std::vector<int> h;


        void moveForward(std::vector<int>& h, ArrayType& particleArray, int u, int n, int m,
            const typename ArrayType::ParticleType& p) {
            int v = h[n + 1] - 1;
            while (m != n) {
                v = h[n + 1] - 1;
                particleArray.changeParticle(particleArray[v], u);
                h[n + 1] = v;
                u = v;
                n = n + 1;
            }
            particleArray.changeParticle(p, v);
        }

        void moveBack(std::vector<int>& h, ArrayType& particleArray, int u, int n, int m,
            const typename ArrayType::ParticleType& p) {
            int v = h[n];
            while (m != n) {
                v = h[n];
                particleArray.changeParticle(particleArray[v], u);
                h[n] = v + 1;
                u = v;
                n = n - 1;
            }
            particleArray.changeParticle(p, v);
        }

        void sortDiap(std::vector<int>& h, ArrayType& particleArray, int n, int u, int s) {
            while (u != s) {
                typename ArrayType::ParticleType p = typename ArrayType::ParticleType(particleArray[u]);
                int m = linIndex(this->grid->getBaseIndex(p.getPosition()));
                if (m == n) {
                    u += 1;
                } 
                else if (m < n) {
                    moveBack(h, particleArray, u, n, m, p);
                    u += 1;
                }
                else {
                    moveForward(h, particleArray, u, n, m, p);
                    s -= 1;
                }
            }
        }
        
        inline int linIndex(const Int3& idx) {
            return (idx.x * this->grid->numCells.y * this->grid->numCells.z +
                    idx.y * this->grid->numCells.z + idx.z);
        }
        inline Int3 cubeIdx(int idx) {
            int x = idx / (this->grid->numCells.y * this->grid->numCells.z);
            idx -= x * (this->grid->numCells.y * this->grid->numCells.z);
            int y = idx / this->grid->numCells.z;
            idx -= y * this->grid->numCells.z;
            return Int3(x, y, idx);
        }
    };
}
