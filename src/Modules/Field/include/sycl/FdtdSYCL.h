#pragma once
#include "Constants.h"
#include "FieldSolver.h"
#include "Grid.h"
#include "Vectors.h"

#include <algorithm>

#include <sycl/CL/sycl.hpp>

#include "../../../../Core/include/sycl/DeviceSYCL.h"
#include "PmlFdtdSYCL.h"

namespace pfc {
    namespace sycl_pfc {
        class FDTD : public RealFieldSolver<YeeGridType>
        {

        public:
            FDTD(YeeGrid* grid, FP dt);

            void updateFields();

            void setPML(int sizePMLx, int sizePMLy, int sizePMLz);
            void setFieldGenerator(FieldGeneratorYee* _generator);

            void updateHalfB();
            void updateE();

            void setTimeStep(FP dt);

            FP getCourantCondition() const {
                double tmp = sqrt(1.0 / (grid->steps.x * grid->steps.x) +
                    1.0 / (grid->steps.y * grid->steps.y) +
                    1.0 / (grid->steps.z * grid->steps.z));
                return 1.0 / (constants::c * tmp);
            }

            bool ifCourantConditionSatisfied(FP dt) const {
                return dt < getCourantCondition();
            }

            void changeDevice(Devices d) {
                device = d;
            }

        private:

            void updateHalfB3D();
            void updateHalfB2D();
            void updateHalfB1D();
            void updateE3D();
            void updateE2D();
            void updateE1D();

            void updateHalfB3D_block();
            void updateE3D_block();


            FP3 anisotropyCoeff;
            void setAnisotropy(const FP frequency, int axis);

            Devices device;

            sycl::range<3> sizeGrid;
            sycl::buffer<FP, 3> bufBx, bufBy, bufBz;
            sycl::buffer<FP, 3> bufEx, bufEy, bufEz;
            sycl::buffer<FP, 3> bufJx, bufJy, bufJz;

            const Int3 block;

        };

        inline FDTD::FDTD(YeeGrid* grid, FP dt) :
            RealFieldSolver(grid, dt, 0.0, 0.5 * dt, 0.0), block(4, 4, 64),
            sizeGrid(grid->sizeStorage.x, grid->sizeStorage.y, grid->sizeStorage.z),
            bufBx(grid->Bx.getData(), sizeGrid),
            bufBy(grid->By.getData(), sizeGrid),
            bufBz(grid->Bz.getData(), sizeGrid),
            bufEx(grid->Ex.getData(), sizeGrid),
            bufEy(grid->Ey.getData(), sizeGrid),
            bufEz(grid->Ez.getData(), sizeGrid),
            bufJx(grid->Jx.getData(), sizeGrid),
            bufJy(grid->Jy.getData(), sizeGrid),
            bufJz(grid->Jz.getData(), sizeGrid)
        {
            if (!ifCourantConditionSatisfied(dt)) {
                std::cout
                    << "WARNING: FDTD Courant condition is not satisfied. Another time step was setted up"
                    << std::endl;
                this->dt = getCourantCondition() * 0.5;
            }
            updateDims();
            pml.reset(new Pml<GridTypes::YeeGridType>(this, Int3(0, 0, 0)));
            generator.reset(new ReflectFieldGeneratorYee(this));
            updateInternalDims();
            anisotropyCoeff = FP3(1, 1, 1);
        }

        inline void FDTD::setPML(int sizePMLx, int sizePMLy, int sizePMLz)
        {
            pml.reset(new PmlFdtd(this, Int3(sizePMLx, sizePMLy, sizePMLz)));
            updateInternalDims();
        }

        inline void FDTD::setTimeStep(FP dt)
        {
            if (ifCourantConditionSatisfied(dt)) {
                this->dt = dt;
                this->timeShiftB = 0.5 * dt;
                if (pml->sizePML == Int3(0, 0, 0))
                    pml.reset(new Pml<GridTypes::YeeGridType>(this, Int3(0, 0, 0)));
                else pml.reset(new PmlFdtd(this, pml->sizePML));
                generator.reset(generator->createInstance(this));
            }
            else {
                std::cout
                    << "WARNING: FDTD Courant condition is not satisfied. Time step was not changed"
                    << std::endl;
            }
        }

        inline void FDTD::setFieldGenerator(FieldGeneratorYee* _generator)
        {
            generator.reset(_generator->createInstance(this));
        }

        inline void FDTD::setAnisotropy(FP frequency, int axis)
        {
            // We introduce artificial anisotropy, through one axis.
            // For this we upgrade Maxwell equations by coefficients,
            // which computes from major signal frequency. 
            // See more in Juntunen,Tsiboukis - Reduction of Numerical Dispersion in
            // FDTD Method Through Artificial Anisotropy.

            FP3 steps = grid->steps;
            FP WP = constants::pi * 2.0 * constants::c / frequency;
            FP R = WP / steps.norm();
            const FP q = 0.99;  // q - stability coefficient, 0 <= q <= 1
            FP Amax = constants::pi / (3 * R * asin(asin(constants::pi / (R * sqrt(3.0))) / sqrt(3.0)));
            FP Q = Amax - 1;
            FP c1 = 1 - Q / 2;
            int axis0 = axis;
            int axis1 = (axis + 1) % 3;
            int axis2 = (axis + 2) % 3;
            // equivalents of the variables
            // Z1 == Zy, Zz == Zz
            // Zy,Zz - designation from article
            FP Z1 = steps[axis0] / steps[axis1];
            FP Z2 = steps[axis0] / steps[axis2];
            // equivalents of the variables
            // CoeffA == K1, CoeffB == K2, a1 == a, a2 == b
            // K1, K2, a, b - designation from article
            FP CoeffA = constants::pi / (R * sqrt(1 + 1 / (Z1 * Z1) + 1 / (Z2 * Z2)));
            FP a1 = sin(CoeffA / c1) * sin(CoeffA / c1)
                / (Z1 * Z1 * sin(CoeffA / (c1 * Z1)) * sin(CoeffA / (c1 * Z1)));
            FP a2 = sin(CoeffA / c1) * sin(CoeffA / c1)
                / (Z2 * Z2 * sin(CoeffA / (c1 * Z2)) * sin(CoeffA / (c1 * Z2)));
            FP CoeffB = sqrt(1 + a1 * Z1 * Z1 + a2 * Z2 * Z2);
            anisotropyCoeff[axis0] = CoeffB / (CoeffA * q * sqrt(a1 * a2))
                * asin(q * sin(CoeffA / c1) / CoeffB);
            anisotropyCoeff[axis1] = a1 * anisotropyCoeff[axis0];
            anisotropyCoeff[axis2] = a2 * anisotropyCoeff[axis0];
        }

        inline void FDTD::updateFields()
        {
            updateHalfB();
            pml->updateB();
            generator->generateB();
            updateE();
            pml->updateE();
            generator->generateE();
            updateHalfB();
            globalTime += dt;
        }

        // Update grid values of magnetic field in FDTD.
        inline void FDTD::updateHalfB()
        {
            if (grid->dimensionality == 3)
                updateHalfB3D_block();
            else if (grid->dimensionality == 2)
                updateHalfB2D();
            else if (grid->dimensionality == 1)
                updateHalfB1D();
        }

        inline void FDTD::updateHalfB3D()
        {
            updateBAreaBegin = Int3(1, 1, 1);
            updateBAreaEnd = grid->numCells - Int3(1, 1, 1);
            for (int d = 0; d < 3; ++d)
            {
                internalBAreaBegin[d] = std::max(updateBAreaBegin[d], pml->leftDims[d]);
                internalBAreaEnd[d] = std::min(updateBAreaEnd[d],
                    grid->numCells[d] - pml->rightDims[d]);
            }

            const FP cdt = constants::c * dt * (FP)0.5;
            const FP coeffXY = cdt / (grid->steps.x * anisotropyCoeff.y);
            const FP coeffXZ = cdt / (grid->steps.x * anisotropyCoeff.z);
            const FP coeffYX = cdt / (grid->steps.y * anisotropyCoeff.x);
            const FP coeffYZ = cdt / (grid->steps.y * anisotropyCoeff.z);
            const FP coeffZX = cdt / (grid->steps.z * anisotropyCoeff.x);
            const FP coeffZY = cdt / (grid->steps.z * anisotropyCoeff.y);

            // In central area use b(i, j, k) += c * dt * -rot(e(i, j, k)), which is:
            // b.x(i, j, k) += c * dt * ((e.y(i, j, k) - e.y(i, j, k-1)) / eps_z * dz -
            //     (e.z(i, j, k) - e.z(i, j-1, k)) / eps_y * dy),
            // b.y(i, j, k) += c * dt * ((e.z(i, j, k) - e.z(i-1, j, k)) / eps_x * dx -
            //     (e.x(i, j, k) - e.x(i, j, k-1)) / eps_z * dz),
            // b.z(i, j, k) += c * dt * ((e.x(i, j, k) - e.x(i, j-1, k)) / eps_y * dy -
            //     (e.y(i, j, k) - e.y(i-1, j, k)) / eps_x * dx),
            try {
                
                sycl::id<3> start(internalBAreaBegin.x, internalBAreaBegin.y, internalBAreaBegin.z);
                sycl::range<3> end(internalBAreaEnd.x, internalBAreaEnd.y, internalBAreaEnd.z);

                auto kernel = [&](sycl::handler& h) {

                    auto Ex = bufEx.get_access<sycl::access::mode::read>(h);
                    auto Ey = bufEy.get_access<sycl::access::mode::read>(h);
                    auto Ez = bufEz.get_access<sycl::access::mode::read>(h);

                    auto Bx = bufBx.get_access<sycl::access::mode::read_write>(h);
                    auto By = bufBy.get_access<sycl::access::mode::read_write>(h);
                    auto Bz = bufBz.get_access<sycl::access::mode::read_write>(h);
                    
                    h.parallel_for(end, start, [=](sycl::item<3> ind) {
                        size_t i = ind[0], j = ind[1], k = ind[2];
                        Bx[i][j][k] += coeffZX * (Ey[i][j][k] - Ey[i][j][k - 1]) -
                            coeffYX * (Ez[i][j][k] - Ez[i][j - 1][k]);
                        By[i][j][k] += coeffXY * (Ez[i][j][k] - Ez[i - 1][j][k]) -
                            coeffZY * (Ex[i][j][k] - Ex[i][j][k - 1]);
                        Bz[i][j][k] += coeffYZ * (Ex[i][j][k] - Ex[i][j - 1][k]) -
                            coeffXZ * (Ey[i][j][k] - Ey[i - 1][j][k]);
                    });
                };

                if (device == Devices::CPU)
                    node.cpu_device.submit(kernel).wait_and_throw();
                if (device == Devices::GPU)
                    node.gpu_device.submit(kernel).wait_and_throw();

                {//update data
                    bufBx.get_access<sycl::access::mode::read_write>();
                    bufBy.get_access<sycl::access::mode::read_write>();
                    bufBz.get_access<sycl::access::mode::read_write>();
                }
            }
            catch (sycl::exception& e)
            {
                std::cout << e.what() << std::endl;
            }
            catch (std::exception& e)
            {
                std::cout << e.what() << std::endl;
            }
        }

        inline void FDTD::updateHalfB3D_block()
        {
            updateBAreaBegin = Int3(1, 1, 1);
            updateBAreaEnd = grid->numCells - Int3(1, 1, 1);
            for (int d = 0; d < 3; ++d)
            {
                internalBAreaBegin[d] = std::max(updateBAreaBegin[d], pml->leftDims[d]);
                internalBAreaEnd[d] = std::min(updateBAreaEnd[d],
                    grid->numCells[d] - pml->rightDims[d]);
            }

            const FP cdt = constants::c * dt * (FP)0.5;
            const FP coeffXY = cdt / (grid->steps.x * anisotropyCoeff.y);
            const FP coeffXZ = cdt / (grid->steps.x * anisotropyCoeff.z);
            const FP coeffYX = cdt / (grid->steps.y * anisotropyCoeff.x);
            const FP coeffYZ = cdt / (grid->steps.y * anisotropyCoeff.z);
            const FP coeffZX = cdt / (grid->steps.z * anisotropyCoeff.x);
            const FP coeffZY = cdt / (grid->steps.z * anisotropyCoeff.y);

            // In central area use b(i, j, k) += c * dt * -rot(e(i, j, k)), which is:
            // b.x(i, j, k) += c * dt * ((e.y(i, j, k) - e.y(i, j, k-1)) / eps_z * dz -
            //     (e.z(i, j, k) - e.z(i, j-1, k)) / eps_y * dy),
            // b.y(i, j, k) += c * dt * ((e.z(i, j, k) - e.z(i-1, j, k)) / eps_x * dx -
            //     (e.x(i, j, k) - e.x(i, j, k-1)) / eps_z * dz),
            // b.z(i, j, k) += c * dt * ((e.x(i, j, k) - e.x(i, j-1, k)) / eps_y * dy -
            //     (e.y(i, j, k) - e.y(i-1, j, k)) / eps_x * dx),
            try {
                
                sycl::id<3> start(internalBAreaBegin.x, internalBAreaBegin.y, internalBAreaBegin.z);
                sycl::id<3> end(internalBAreaEnd.x, internalBAreaEnd.y, internalBAreaEnd.z);
                sycl::range<3> data_size(block.x + 1, block.y + 1, block.z + 1);
                sycl::range<3> block_size(block.x, block.y, block.z);
                sycl::range<3> block_count((internalBAreaEnd.x - internalBAreaBegin.x + block.x - 1) / block.x,
                                           (internalBAreaEnd.y - internalBAreaBegin.y + block.y - 1) / block.y,
                                           (internalBAreaEnd.z - internalBAreaBegin.z + block.z - 1) / block.z);


                auto kernel = [&](sycl::handler& h) {

                    auto Ex = bufEx.get_access<sycl::access::mode::read>(h);
                    auto Ey = bufEy.get_access<sycl::access::mode::read>(h);
                    auto Ez = bufEz.get_access<sycl::access::mode::read>(h);

                    auto Bx = bufBx.get_access<sycl::access::mode::read_write>(h);
                    auto By = bufBy.get_access<sycl::access::mode::read_write>(h);
                    auto Bz = bufBz.get_access<sycl::access::mode::read_write>(h);

                    using LocalAccessor3d = sycl::accessor<FP, 3, 
                          sycl::access::mode::read_write, sycl::access::target::local>;
                    
                    LocalAccessor3d Ex_local(data_size, h);
                    LocalAccessor3d Ey_local(data_size, h);
                    LocalAccessor3d Ez_local(data_size, h);

                    h.parallel_for_work_group(block_count, [=](sycl::group<3> grp) {
                        auto grpInd = grp.get_id();
                        sycl::id start_block(grpInd[0] * block_size[0], 
                                             grpInd[1] * block_size[1], grpInd[2] * block_size[2]);
                        start_block += start;
 
                        Int3 n_block;
                        n_block[0] = start_block[0] + block_size[0] < end[0] ? 
                               block_size[0] : end[0] - start_block[0];
                        n_block[1] = start_block[1] + block_size[1] < end[1] ? 
                               block_size[1] : end[1] - start_block[1];
                        n_block[2] = start_block[2] + block_size[2] < end[2] ? 
                               block_size[2] : end[2] - start_block[2];
                        Int3 n_data = n_block + Int3(1, 1, 1);

                        grp.parallel_for_work_item(sycl::range(n_data.x, n_data.y, n_data.z),
                           [=](sycl::h_item<3> myItem){
                             auto it_local = myItem.get_local_id();
                             auto it = start_block + myItem.get_local_id() - sycl::id(1, 1, 1);
                             //if(it[0] < end[0] && it[1] < end[1] && it[2] < end[2])
                             {
                                 Ex_local[it_local] = Ex[it];
                                 Ey_local[it_local] = Ey[it];
                                 Ez_local[it_local] = Ez[it];
                             }
                        });

                        grp.parallel_for_work_item(sycl::range(n_block.x, n_block.y, n_block.z),
                           [=](sycl::h_item<3> myItem){
                             auto it_local = myItem.get_local_id();
                             auto it = start_block + myItem.get_local_id();
                             size_t i = it_local[0] + 1, j = it_local[1] + 1, k = it_local[2] + 1;

                             Bx[it] += coeffZX * (Ey_local[i][j][k] - Ey_local[i][j][k - 1]) -
                                 coeffYX * (Ez_local[i][j][k] - Ez_local[i][j - 1][k]);
                             By[it] += coeffXY * (Ez_local[i][j][k] - Ez_local[i - 1][j][k]) -
                                 coeffZY * (Ex_local[i][j][k] - Ex_local[i][j][k - 1]);
                             Bz[it] += coeffYZ * (Ex_local[i][j][k] - Ex_local[i][j - 1][k]) -
                                 coeffXZ * (Ey_local[i][j][k] - Ey_local[i - 1][j][k]);

                        });

                    });

                    /*h.parallel_for(end, start, [=](sycl::item<3> ind) {
                        size_t i = ind[0], j = ind[1], k = ind[2];
                        Bx[i][j][k] += coeffZX * (Ey[i][j][k] - Ey[i][j][k - 1]) -
                            coeffYX * (Ez[i][j][k] - Ez[i][j - 1][k]);
                        By[i][j][k] += coeffXY * (Ez[i][j][k] - Ez[i - 1][j][k]) -
                            coeffZY * (Ex[i][j][k] - Ex[i][j][k - 1]);
                        Bz[i][j][k] += coeffYZ * (Ex[i][j][k] - Ex[i][j - 1][k]) -
                            coeffXZ * (Ey[i][j][k] - Ey[i - 1][j][k]);
                    });*/
                };

                if (device == Devices::CPU)
                    node.cpu_device.submit(kernel).wait_and_throw();
                if (device == Devices::GPU)
                    node.gpu_device.submit(kernel).wait_and_throw();

                {//update data
                    bufBx.get_access<sycl::access::mode::read_write>();
                    bufBy.get_access<sycl::access::mode::read_write>();
                    bufBz.get_access<sycl::access::mode::read_write>();
                }
            }
            catch (sycl::exception& e)
            {
                std::cout << e.what() << std::endl;
            }
            catch (std::exception& e)
            {
                std::cout << e.what() << std::endl;
            }
        }

        inline void FDTD::updateHalfB2D()
        {
            updateBAreaBegin = Int3(1, 1, 0);
            updateBAreaEnd = grid->numCells - Int3(1, 1, 0);
            for (int d = 0; d < 2; ++d)
            {
                internalBAreaBegin[d] = std::max(updateBAreaBegin[d], pml->leftDims[d]);
                internalBAreaEnd[d] = std::min(updateBAreaEnd[d],
                    grid->numCells[d] - pml->rightDims[d]);
            }

            const FP cdt = constants::c * dt * (FP)0.5;
            const FP coeffXY = cdt / (grid->steps.x * anisotropyCoeff.y);
            const FP coeffXZ = cdt / (grid->steps.x * anisotropyCoeff.z);
            const FP coeffYX = cdt / (grid->steps.y * anisotropyCoeff.x);
            const FP coeffYZ = cdt / (grid->steps.y * anisotropyCoeff.z);

            // In central area use b(i, j, k) += c * dt * -rot(e(i, j, k)), which is:
            // b.x(i, j, k) += c * dt * ((e.y(i, j, k) - e.y(i, j, k-1)) / eps_z * dz -
            //     (e.z(i, j, k) - e.z(i, j-1, k)) / eps_y * dy),
            // b.y(i, j, k) += c * dt * ((e.z(i, j, k) - e.z(i-1, j, k)) / eps_x * dx -
            //     (e.x(i, j, k) - e.x(i, j, k-1)) / eps_z * dz),
            // b.z(i, j, k) += c * dt * ((e.x(i, j, k) - e.x(i, j-1, k)) / eps_y * dy -
            //     (e.y(i, j, k) - e.y(i-1, j, k)) / eps_x * dx),
            try {

                sycl::id<2> start(internalBAreaBegin.x, internalBAreaBegin.y);
                sycl::range<2> end(internalBAreaEnd.x, internalBAreaEnd.y);

                auto kernel = [&](sycl::handler& h) {

                    auto Ex = bufEx.get_access<sycl::access::mode::read>(h);
                    auto Ey = bufEy.get_access<sycl::access::mode::read>(h);
                    auto Ez = bufEz.get_access<sycl::access::mode::read>(h);

                    auto Bx = bufBx.get_access<sycl::access::mode::read_write>(h);
                    auto By = bufBy.get_access<sycl::access::mode::read_write>(h);
                    auto Bz = bufBz.get_access<sycl::access::mode::read_write>(h);

                    h.parallel_for(end, start, [=](sycl::item<2> ind) {
                        size_t i = ind[0], j = ind[1], k = 0;
                        Bx[i][j][k] += -coeffYX * (Ez[i][j][k] - Ez[i][j - 1][k]);
                        By[i][j][k] += coeffXY * (Ez[i][j][k] - Ez[i - 1][j][k]);
                        Bz[i][j][k] += coeffYZ * (Ex[i][j][k] - Ex[i][j - 1][k]) -
                            coeffXZ * (Ey[i][j][k] - Ey[i - 1][j][k]);
                    });
                };

                if (device == Devices::CPU)
                    node.cpu_device.submit(kernel).wait_and_throw();
                if (device == Devices::GPU)
                    node.gpu_device.submit(kernel).wait_and_throw();

                {//update data
                    bufBx.get_access<sycl::access::mode::read_write>();
                    bufBy.get_access<sycl::access::mode::read_write>();
                    bufBz.get_access<sycl::access::mode::read_write>();
                }
            }
            catch (sycl::exception& e)
            {
                std::cout << e.what() << std::endl;
            }
            catch (std::exception& e)
            {
                std::cout << e.what() << std::endl;
            }
        }

        inline void FDTD::updateHalfB1D()
        {
            updateBAreaBegin = Int3(1, 0, 0);
            updateBAreaEnd = grid->numCells - Int3(1, 0, 0);
            for (int d = 0; d < 1; ++d)
            {
                internalBAreaBegin[d] = std::max(updateBAreaBegin[d], pml->leftDims[d]);
                internalBAreaEnd[d] = std::min(updateBAreaEnd[d],
                    grid->numCells[d] - pml->rightDims[d]);
            }

            const FP cdt = constants::c * dt * (FP)0.5;
            const FP coeffXY = cdt / (grid->steps.x * anisotropyCoeff.y);
            const FP coeffXZ = cdt / (grid->steps.x * anisotropyCoeff.z);

            // In central area use b(i, j, k) += c * dt * -rot(e(i, j, k)), which is:
            // b.x(i, j, k) += c * dt * ((e.y(i, j, k) - e.y(i, j, k-1)) / eps_z * dz -
            //     (e.z(i, j, k) - e.z(i, j-1, k)) / eps_y * dy),
            // b.y(i, j, k) += c * dt * ((e.z(i, j, k) - e.z(i-1, j, k)) / eps_x * dx -
            //     (e.x(i, j, k) - e.x(i, j, k-1)) / eps_z * dz),
            // b.z(i, j, k) += c * dt * ((e.x(i, j, k) - e.x(i, j-1, k)) / eps_y * dy -
            //     (e.y(i, j, k) - e.y(i-1, j, k)) / eps_x * dx),
            try {

                sycl::id<1> start(internalBAreaBegin.x);
                sycl::range<1> end(internalBAreaEnd.x);

                auto kernel = [&](sycl::handler& h) {

                    auto Ex = bufEx.get_access<sycl::access::mode::read>(h);
                    auto Ey = bufEy.get_access<sycl::access::mode::read>(h);
                    auto Ez = bufEz.get_access<sycl::access::mode::read>(h);

                    auto By = bufBy.get_access<sycl::access::mode::read_write>(h);
                    auto Bz = bufBz.get_access<sycl::access::mode::read_write>(h);

                    h.parallel_for(end, start, [=](sycl::item<1> ind) {
                        size_t i = ind[0], j = 0, k = 0;
                        By[i][j][k] += coeffXY * (Ez[i][j][k] - Ez[i - 1][j][k]);
                        Bz[i][j][k] += -coeffXZ * (Ey[i][j][k] - Ey[i - 1][j][k]);
                    });
                };

                if (device == Devices::CPU)
                    node.cpu_device.submit(kernel).wait_and_throw();
                if (device == Devices::GPU)
                    node.gpu_device.submit(kernel).wait_and_throw();

                {//update data
                    bufBy.get_access<sycl::access::mode::read_write>();
                    bufBz.get_access<sycl::access::mode::read_write>();
                }

            }
            catch (sycl::exception& e)
            {
                std::cout << e.what() << std::endl;
            }
            catch (std::exception& e)
            {
                std::cout << e.what() << std::endl;
            }
        }

        // Update grid values of electric field in FDTD.
        inline void FDTD::updateE()
        {
            if (grid->dimensionality == 3)
                updateE3D_block();
            else if (grid->dimensionality == 2)
                updateE2D();
            else if (grid->dimensionality == 1)
                updateE1D();
        }

        inline void FDTD::updateE3D()
        {
            updateEAreaBegin = Int3(0, 0, 0);
            updateEAreaEnd = grid->numCells - Int3(1, 1, 1);
            for (int d = 0; d < 3; ++d)
            {
                internalEAreaBegin[d] = std::max(updateEAreaBegin[d], pml->leftDims[d]);
                internalEAreaEnd[d] = std::min(updateEAreaEnd[d],
                    grid->numCells[d] - pml->rightDims[d]);
            }

            const FP coeffCurrent = -(FP)4 * constants::pi * dt;
            const FP cdt = constants::c * dt;
            const FP coeffXY = cdt / (grid->steps.x * anisotropyCoeff.y);
            const FP coeffXZ = cdt / (grid->steps.x * anisotropyCoeff.z);
            const FP coeffYX = cdt / (grid->steps.y * anisotropyCoeff.x);
            const FP coeffYZ = cdt / (grid->steps.y * anisotropyCoeff.z);
            const FP coeffZX = cdt / (grid->steps.z * anisotropyCoeff.x);
            const FP coeffZY = cdt / (grid->steps.z * anisotropyCoeff.y);

            // In internal area use:
            // e.x(i, j, k) += dt * -4pi * j.x(i, j, k) + c * dt * ((b.z(i, j+1, k) -
            //     b.z(i, j, k)) / eps_y * dy - (b.y(i, j, k+1) - b.y(i, j, k)) / eps_z * dz),
            // e.y(i, j, k) += dt * -4pi * j.y(i, j, k) + c * dt * ((b.x(i, j, k+1) -
            //     b.x(i, j, k)) / eps_z * dz - (b.z(i+1, j, k) - b.z(i, j, k)) / eps_x * dx),
            // e.z(i, j, k) += dt * -4pi * j.z(i, j, k) + c * dt * ((b.y(i+1, j, k) -
            //     b.y(i, j, k)) / eps_x * dx - (b.x(i, j+1, k) - b.x(i, j, k)) / eps_y * dy),
            try {
                
                sycl::id<3> start(internalEAreaBegin.x, internalEAreaBegin.y, internalEAreaBegin.z);
                sycl::range<3> end(internalEAreaEnd.x, internalEAreaEnd.y, internalEAreaEnd.z);

                auto kernel = [&](sycl::handler& h) {

                    auto Bx = bufBx.get_access<sycl::access::mode::read>(h);
                    auto By = bufBy.get_access<sycl::access::mode::read>(h);
                    auto Bz = bufBz.get_access<sycl::access::mode::read>(h);
                    auto Jx = bufJx.get_access<sycl::access::mode::read>(h);
                    auto Jy = bufJy.get_access<sycl::access::mode::read>(h);
                    auto Jz = bufJz.get_access<sycl::access::mode::read>(h);

                    auto Ex = bufEx.get_access<sycl::access::mode::read_write>(h);
                    auto Ey = bufEy.get_access<sycl::access::mode::read_write>(h);
                    auto Ez = bufEz.get_access<sycl::access::mode::read_write>(h);

                    h.parallel_for(end, start, [=](sycl::item<3> ind) {
                        size_t i = ind[0], j = ind[1], k = ind[2];
                        Ex[i][j][k] += coeffCurrent * Jx[i][j][k] +
                            coeffYX * (Bz[i][j + 1][k] - Bz[i][j][k]) -
                            coeffZX * (By[i][j][k + 1] - By[i][j][k]);
                        Ey[i][j][k] += coeffCurrent * Jy[i][j][k] +
                            coeffZY * (Bx[i][j][k + 1] - Bx[i][j][k]) -
                            coeffXY * (Bz[i + 1][j][k] - Bz[i][j][k]);
                        Ez[i][j][k] += coeffCurrent * Jz[i][j][k] +
                            coeffXZ * (By[i + 1][j][k] - By[i][j][k]) -
                            coeffYZ * (Bx[i][j + 1][k] - Bx[i][j][k]);
                    });
                };

                // Process edge values

                auto kernel_boundX = [&](sycl::handler& h) {

                    auto By = bufBy.get_access<sycl::access::mode::read>(h);
                    auto Bz = bufBz.get_access<sycl::access::mode::read>(h);
                    auto Jx = bufJx.get_access<sycl::access::mode::read>(h);
                    
                    auto Ex = bufEx.get_access<sycl::access::mode::read_write>(h);
                    
                    size_t i = updateEAreaEnd.x;
                    h.parallel_for(sycl::range(internalEAreaEnd.y, internalEAreaEnd.z),
                        sycl::id<2>(internalEAreaBegin.y, internalEAreaBegin.z),
                        [=](sycl::item<2> ind) {
                        size_t j = ind[0], k = ind[1];
                        Ex[i][j][k] += coeffCurrent * Jx[i][j][k] +
                            coeffYX * (Bz[i][j + 1][k] - Bz[i][j][k]) -
                            coeffZX * (By[i][j][k + 1] - By[i][j][k]);
                    });
                };
                auto kernel_boundY = [&](sycl::handler& h) {

                    auto Bx = bufBx.get_access<sycl::access::mode::read>(h);
                    auto Bz = bufBz.get_access<sycl::access::mode::read>(h);
                    auto Jy = bufJy.get_access<sycl::access::mode::read>(h);
                    
                    auto Ey = bufEy.get_access<sycl::access::mode::read_write>(h);
                    
                    size_t j = updateEAreaEnd.y;
                    h.parallel_for(sycl::range(internalEAreaEnd.x, internalEAreaEnd.z),
                        sycl::id<2>(internalEAreaBegin.x, internalEAreaBegin.z),
                        [=](sycl::item<2> ind) {
                        size_t i = ind[0], k = ind[1];
                        Ey[i][j][k] += coeffCurrent * Jy[i][j][k] +
                            coeffZY * (Bx[i][j][k + 1] - Bx[i][j][k]) -
                            coeffXY * (Bz[i + 1][j][k] - Bz[i][j][k]);
                    });
                };
                auto kernel_boundZ = [&](sycl::handler& h) {

                    auto Bx = bufBx.get_access<sycl::access::mode::read>(h);
                    auto By = bufBy.get_access<sycl::access::mode::read>(h);
                    auto Jz = bufJz.get_access<sycl::access::mode::read>(h);

                    auto Ez = bufEz.get_access<sycl::access::mode::read_write>(h);

                    size_t k = updateEAreaEnd.z;
                    h.parallel_for(sycl::range(internalEAreaEnd.x, internalEAreaEnd.y),
                        sycl::id<2>(internalEAreaBegin.x, internalEAreaBegin.y),
                        [=](sycl::item<2> ind) {
                        size_t i = ind[0], j = ind[1];
                        Ez[i][j][k] += coeffCurrent * Jz[i][j][k] +
                            coeffXZ * (By[i + 1][j][k] - By[i][j][k]) -
                            coeffYZ * (Bx[i][j + 1][k] - Bx[i][j][k]);
                    });
                };

                if (device == Devices::CPU)
                {
                    node.cpu_device.submit(kernel).wait_and_throw();
                    if (updateEAreaEnd.x == grid->numCells.x - 1)
                        node.cpu_device.submit(kernel_boundX).wait_and_throw();
                    if (updateEAreaEnd.y == grid->numCells.y - 1)
                        node.cpu_device.submit(kernel_boundY).wait_and_throw();
                    if (updateEAreaEnd.z == grid->numCells.z - 1)
                        node.cpu_device.submit(kernel_boundZ).wait_and_throw();
                }
                if (device == Devices::GPU)
                {
                    node.gpu_device.submit(kernel).wait_and_throw();
                    if (updateEAreaEnd.x == grid->numCells.x - 1)
                        node.gpu_device.submit(kernel_boundX).wait_and_throw();
                    if (updateEAreaEnd.y == grid->numCells.y - 1)
                        node.gpu_device.submit(kernel_boundY).wait_and_throw();
                    if (updateEAreaEnd.z == grid->numCells.z - 1)
                        node.gpu_device.submit(kernel_boundZ).wait_and_throw();
                }

                {//update data
                    bufEx.get_access<sycl::access::mode::read_write>();
                    bufEy.get_access<sycl::access::mode::read_write>();
                    bufEz.get_access<sycl::access::mode::read_write>();
                }
            }
            catch (sycl::exception& e)
            {
                std::cout << e.what() << std::endl;
            }
            catch (std::exception& e)
            {
                std::cout << e.what() << std::endl;
            }
        }

        inline void FDTD::updateE3D_block()
        {
            updateEAreaBegin = Int3(0, 0, 0);
            updateEAreaEnd = grid->numCells - Int3(1, 1, 1);
            for (int d = 0; d < 3; ++d)
            {
                internalEAreaBegin[d] = std::max(updateEAreaBegin[d], pml->leftDims[d]);
                internalEAreaEnd[d] = std::min(updateEAreaEnd[d],
                    grid->numCells[d] - pml->rightDims[d]);
            }

            const FP coeffCurrent = -(FP)4 * constants::pi * dt;
            const FP cdt = constants::c * dt;
            const FP coeffXY = cdt / (grid->steps.x * anisotropyCoeff.y);
            const FP coeffXZ = cdt / (grid->steps.x * anisotropyCoeff.z);
            const FP coeffYX = cdt / (grid->steps.y * anisotropyCoeff.x);
            const FP coeffYZ = cdt / (grid->steps.y * anisotropyCoeff.z);
            const FP coeffZX = cdt / (grid->steps.z * anisotropyCoeff.x);
            const FP coeffZY = cdt / (grid->steps.z * anisotropyCoeff.y);

            // In internal area use:
            // e.x(i, j, k) += dt * -4pi * j.x(i, j, k) + c * dt * ((b.z(i, j+1, k) -
            //     b.z(i, j, k)) / eps_y * dy - (b.y(i, j, k+1) - b.y(i, j, k)) / eps_z * dz),
            // e.y(i, j, k) += dt * -4pi * j.y(i, j, k) + c * dt * ((b.x(i, j, k+1) -
            //     b.x(i, j, k)) / eps_z * dz - (b.z(i+1, j, k) - b.z(i, j, k)) / eps_x * dx),
            // e.z(i, j, k) += dt * -4pi * j.z(i, j, k) + c * dt * ((b.y(i+1, j, k) -
            //     b.y(i, j, k)) / eps_x * dx - (b.x(i, j+1, k) - b.x(i, j, k)) / eps_y * dy),
            try {
                
                sycl::id<3> start(internalEAreaBegin.x, internalEAreaBegin.y, internalEAreaBegin.z);
                sycl::id<3> end(internalEAreaEnd.x, internalEAreaEnd.y, internalEAreaEnd.z);
                sycl::range<3> data_size(block.x + 1, block.y + 1, block.z + 1);
                sycl::range<3> block_size(block.x, block.y, block.z);
                sycl::range<3> block_count((internalEAreaEnd.x - internalEAreaBegin.x + block.x - 1) / block.x,
                                           (internalEAreaEnd.y - internalEAreaBegin.y + block.y - 1) / block.y,
                                           (internalEAreaEnd.z - internalEAreaBegin.z + block.z - 1) / block.z);

                auto kernel = [&](sycl::handler& h) {

                    auto Bx = bufBx.get_access<sycl::access::mode::read>(h);
                    auto By = bufBy.get_access<sycl::access::mode::read>(h);
                    auto Bz = bufBz.get_access<sycl::access::mode::read>(h);
                    auto Jx = bufJx.get_access<sycl::access::mode::read>(h);
                    auto Jy = bufJy.get_access<sycl::access::mode::read>(h);
                    auto Jz = bufJz.get_access<sycl::access::mode::read>(h);

                    auto Ex = bufEx.get_access<sycl::access::mode::read_write>(h);
                    auto Ey = bufEy.get_access<sycl::access::mode::read_write>(h);
                    auto Ez = bufEz.get_access<sycl::access::mode::read_write>(h);

                    using LocalAccessor3d = sycl::accessor<FP, 3, 
                          sycl::access::mode::read_write, sycl::access::target::local>;
                    
                    LocalAccessor3d Bx_local(data_size, h);
                    LocalAccessor3d By_local(data_size, h);
                    LocalAccessor3d Bz_local(data_size, h);

                    h.parallel_for_work_group(block_count, [=](sycl::group<3> grp) {
                        auto grpInd = grp.get_id();
                        sycl::id start_block(grpInd[0] * block_size[0], 
                                             grpInd[1] * block_size[1], grpInd[2] * block_size[2]);
                        start_block += start;

                        grp.parallel_for_work_item(data_size, [=](sycl::h_item<3> myItem){
                             auto it_local = myItem.get_local_id();
                             auto it = start_block + myItem.get_local_id();
                             if(it[0] < end[0] && it[1] < end[1] && it[2] < end[2]){
                                 Bx_local[it_local] = Bx[it];
                                 By_local[it_local] = By[it];
                                 Bz_local[it_local] = Bz[it];
                             }
                        });

                        grp.parallel_for_work_item(block_size, [=](sycl::h_item<3> myItem){
                             auto it_local = myItem.get_local_id();
                             auto it = start_block + myItem.get_local_id();
                             size_t i = it_local[0], j = it_local[1], k = it_local[2];

                             Ex[it] += coeffCurrent * Jx[it] +
                                       coeffYX * (Bz_local[i][j + 1][k] - Bz_local[i][j][k]) -
                                       coeffZX * (By_local[i][j][k + 1] - By_local[i][j][k]);
                             Ey[it] += coeffCurrent * Jy[it] +
                                       coeffZY * (Bx_local[i][j][k + 1] - Bx_local[i][j][k]) -
                                       coeffXY * (Bz_local[i + 1][j][k] - Bz_local[i][j][k]);
                             Ez[it] += coeffCurrent * Jz[it] +
                                       coeffXZ * (By_local[i + 1][j][k] - By_local[i][j][k]) -
                                       coeffYZ * (Bx_local[i][j + 1][k] - Bx_local[i][j][k]);

                        });

                    });

                    /*h.parallel_for(end, start, [=](sycl::item<3> ind) {
                        size_t i = ind[0], j = ind[1], k = ind[2];
                        Ex[i][j][k] += coeffCurrent * Jx[i][j][k] +
                            coeffYX * (Bz[i][j + 1][k] - Bz[i][j][k]) -
                            coeffZX * (By[i][j][k + 1] - By[i][j][k]);
                        Ey[i][j][k] += coeffCurrent * Jy[i][j][k] +
                            coeffZY * (Bx[i][j][k + 1] - Bx[i][j][k]) -
                            coeffXY * (Bz[i + 1][j][k] - Bz[i][j][k]);
                        Ez[i][j][k] += coeffCurrent * Jz[i][j][k] +
                            coeffXZ * (By[i + 1][j][k] - By[i][j][k]) -
                            coeffYZ * (Bx[i][j + 1][k] - Bx[i][j][k]);
                    });*/
                };

                // Process edge values

                auto kernel_boundX = [&](sycl::handler& h) {

                    auto By = bufBy.get_access<sycl::access::mode::read>(h);
                    auto Bz = bufBz.get_access<sycl::access::mode::read>(h);
                    auto Jx = bufJx.get_access<sycl::access::mode::read>(h);
                    
                    auto Ex = bufEx.get_access<sycl::access::mode::read_write>(h);
                    
                    size_t i = updateEAreaEnd.x;
                    h.parallel_for(sycl::range(internalEAreaEnd.y, internalEAreaEnd.z),
                        sycl::id<2>(internalEAreaBegin.y, internalEAreaBegin.z),
                        [=](sycl::item<2> ind) {
                        size_t j = ind[0], k = ind[1];
                        Ex[i][j][k] += coeffCurrent * Jx[i][j][k] +
                            coeffYX * (Bz[i][j + 1][k] - Bz[i][j][k]) -
                            coeffZX * (By[i][j][k + 1] - By[i][j][k]);
                    });
                };
                auto kernel_boundY = [&](sycl::handler& h) {

                    auto Bx = bufBx.get_access<sycl::access::mode::read>(h);
                    auto Bz = bufBz.get_access<sycl::access::mode::read>(h);
                    auto Jy = bufJy.get_access<sycl::access::mode::read>(h);
                    
                    auto Ey = bufEy.get_access<sycl::access::mode::read_write>(h);
                    
                    size_t j = updateEAreaEnd.y;
                    h.parallel_for(sycl::range(internalEAreaEnd.x, internalEAreaEnd.z),
                        sycl::id<2>(internalEAreaBegin.x, internalEAreaBegin.z),
                        [=](sycl::item<2> ind) {
                        size_t i = ind[0], k = ind[1];
                        Ey[i][j][k] += coeffCurrent * Jy[i][j][k] +
                            coeffZY * (Bx[i][j][k + 1] - Bx[i][j][k]) -
                            coeffXY * (Bz[i + 1][j][k] - Bz[i][j][k]);
                    });
                };
                auto kernel_boundZ = [&](sycl::handler& h) {

                    auto Bx = bufBx.get_access<sycl::access::mode::read>(h);
                    auto By = bufBy.get_access<sycl::access::mode::read>(h);
                    auto Jz = bufJz.get_access<sycl::access::mode::read>(h);

                    auto Ez = bufEz.get_access<sycl::access::mode::read_write>(h);

                    size_t k = updateEAreaEnd.z;
                    h.parallel_for(sycl::range(internalEAreaEnd.x, internalEAreaEnd.y),
                        sycl::id<2>(internalEAreaBegin.x, internalEAreaBegin.y),
                        [=](sycl::item<2> ind) {
                        size_t i = ind[0], j = ind[1];
                        Ez[i][j][k] += coeffCurrent * Jz[i][j][k] +
                            coeffXZ * (By[i + 1][j][k] - By[i][j][k]) -
                            coeffYZ * (Bx[i][j + 1][k] - Bx[i][j][k]);
                    });
                };

                if (device == Devices::CPU)
                {
                    node.cpu_device.submit(kernel).wait_and_throw();
                    if (updateEAreaEnd.x == grid->numCells.x - 1)
                        node.cpu_device.submit(kernel_boundX).wait_and_throw();
                    if (updateEAreaEnd.y == grid->numCells.y - 1)
                        node.cpu_device.submit(kernel_boundY).wait_and_throw();
                    if (updateEAreaEnd.z == grid->numCells.z - 1)
                        node.cpu_device.submit(kernel_boundZ).wait_and_throw();
                }
                if (device == Devices::GPU)
                {
                    node.gpu_device.submit(kernel).wait_and_throw();
                    if (updateEAreaEnd.x == grid->numCells.x - 1)
                        node.gpu_device.submit(kernel_boundX).wait_and_throw();
                    if (updateEAreaEnd.y == grid->numCells.y - 1)
                        node.gpu_device.submit(kernel_boundY).wait_and_throw();
                    if (updateEAreaEnd.z == grid->numCells.z - 1)
                        node.gpu_device.submit(kernel_boundZ).wait_and_throw();
                }

                {//update data
                    bufEx.get_access<sycl::access::mode::read_write>();
                    bufEy.get_access<sycl::access::mode::read_write>();
                    bufEz.get_access<sycl::access::mode::read_write>();
                }
            }
            catch (sycl::exception& e)
            {
                std::cout << e.what() << std::endl;
            }
            catch (std::exception& e)
            {
                std::cout << e.what() << std::endl;
            }
        }

        inline void FDTD::updateE2D()
        {
            updateEAreaBegin = Int3(0, 0, 0);
            updateEAreaEnd = grid->numCells - Int3(1, 1, 0);
            for (int d = 0; d < 2; ++d)
            {
                internalEAreaBegin[d] = std::max(updateEAreaBegin[d], pml->leftDims[d]);
                internalEAreaEnd[d] = std::min(updateEAreaEnd[d],
                    grid->numCells[d] - pml->rightDims[d]);
            }

            const FP coeffCurrent = -(FP)4 * constants::pi * dt;
            const FP cdt = constants::c * dt;
            const FP coeffXY = cdt / (grid->steps.x * anisotropyCoeff.y);
            const FP coeffXZ = cdt / (grid->steps.x * anisotropyCoeff.z);
            const FP coeffYX = cdt / (grid->steps.y * anisotropyCoeff.x);
            const FP coeffYZ = cdt / (grid->steps.y * anisotropyCoeff.z);

            // In internal area use:
            // e.x(i, j, k) += dt * -4pi * j.x(i, j, k) + c * dt * ((b.z(i, j+1, k) -
            //     b.z(i, j, k)) / eps_y * dy - (b.y(i, j, k+1) - b.y(i, j, k)) / eps_z * dz),
            // e.y(i, j, k) += dt * -4pi * j.y(i, j, k) + c * dt * ((b.x(i, j, k+1) -
            //     b.x(i, j, k)) / eps_z * dz - (b.z(i+1, j, k) - b.z(i, j, k)) / eps_x * dx),
            // e.z(i, j, k) += dt * -4pi * j.z(i, j, k) + c * dt * ((b.y(i+1, j, k) -
            //     b.y(i, j, k)) / eps_x * dx - (b.x(i, j+1, k) - b.x(i, j, k)) / eps_y * dy),
            try {

                sycl::id<2> start(internalEAreaBegin.x, internalEAreaBegin.y);
                sycl::range<2> end(internalEAreaEnd.x, internalEAreaEnd.y);

                auto kernel = [&](sycl::handler& h) {

                    auto Bx = bufBx.get_access<sycl::access::mode::read>(h);
                    auto By = bufBy.get_access<sycl::access::mode::read>(h);
                    auto Bz = bufBz.get_access<sycl::access::mode::read>(h);
                    auto Jx = bufJx.get_access<sycl::access::mode::read>(h);
                    auto Jy = bufJy.get_access<sycl::access::mode::read>(h);
                    auto Jz = bufJz.get_access<sycl::access::mode::read>(h);

                    auto Ex = bufEx.get_access<sycl::access::mode::read_write>(h);
                    auto Ey = bufEy.get_access<sycl::access::mode::read_write>(h);
                    auto Ez = bufEz.get_access<sycl::access::mode::read_write>(h);

                    h.parallel_for(end, start, [=](sycl::item<2> ind) {
                        size_t i = ind[0], j = ind[1], k = 0;
                        Ex[i][j][k] += coeffCurrent * Jx[i][j][k] +
                            coeffYX * (Bz[i][j + 1][k] - Bz[i][j][k]);
                        Ey[i][j][k] += coeffCurrent * Jy[i][j][k] -
                            coeffXY * (Bz[i + 1][j][k] - Bz[i][j][k]);
                        Ez[i][j][k] += coeffCurrent * Jz[i][j][k] +
                            coeffXZ * (By[i + 1][j][k] - By[i][j][k]) -
                            coeffYZ * (Bx[i][j + 1][k] - Bx[i][j][k]);
                    });
                };
                auto kernel_boundX = [&](sycl::handler& h) {

                    auto Bz = bufBz.get_access<sycl::access::mode::read>(h);
                    auto Jx = bufJx.get_access<sycl::access::mode::read>(h);
                    
                    auto Ex = bufEx.get_access<sycl::access::mode::read_write>(h);
                    
                    size_t i = updateEAreaEnd.x, k = 0;
                    h.parallel_for(sycl::range(internalEAreaEnd.y),
                        sycl::id<1>(internalEAreaBegin.y),
                        [=](sycl::item<1> ind) {
                        size_t j = ind[0];
                        Ex[i][j][k] += coeffCurrent * Jx[i][j][k] +
                            coeffYX * (Bz[i][j + 1][k] - Bz[i][j][k]);
                    });
                };
                auto kernel_boundY = [&](sycl::handler& h) {

                    auto Bz = bufBz.get_access<sycl::access::mode::read>(h);
                    auto Jy = bufJy.get_access<sycl::access::mode::read>(h);
                    
                    auto Ey = bufEy.get_access<sycl::access::mode::read_write>(h);
                    
                    size_t j = updateEAreaEnd.y, k = 0;
                    h.parallel_for(sycl::range(internalEAreaEnd.x),
                        sycl::id<1>(internalEAreaBegin.x),
                        [=](sycl::item<1> ind) {
                        size_t i = ind[0];
                        Ey[i][j][k] += coeffCurrent * Jy[i][j][k] -
                            coeffXY * (Bz[i + 1][j][k] - Bz[i][j][k]);
                    });
                };

                if (device == Devices::CPU)
                {
                    node.cpu_device.submit(kernel);
                    if (updateEAreaEnd.x == grid->numCells.x - 1)
                        node.cpu_device.submit(kernel_boundX);
                    if (updateEAreaEnd.y == grid->numCells.y - 1)
                        node.cpu_device.submit(kernel_boundY);
                    node.cpu_device.wait_and_throw();
                }
                if (device == Devices::GPU)
                {
                    node.gpu_device.submit(kernel).wait_and_throw();
                    if (updateEAreaEnd.x == grid->numCells.x - 1)
                        node.gpu_device.submit(kernel_boundX).wait_and_throw();
                    if (updateEAreaEnd.y == grid->numCells.y - 1)
                        node.gpu_device.submit(kernel_boundY).wait_and_throw();
                }

                {//update data
                    bufEx.get_access<sycl::access::mode::read_write>();
                    bufEy.get_access<sycl::access::mode::read_write>();
                    bufEz.get_access<sycl::access::mode::read_write>();
                }
            }
            catch (sycl::exception& e)
            {
                std::cout << e.what() << std::endl;
            }
            catch (std::exception& e)
            {
                std::cout << e.what() << std::endl;
            }
        }

        inline void FDTD::updateE1D()
        {
            updateEAreaBegin = Int3(0, 0, 0);
            updateEAreaEnd = grid->numCells - Int3(1, 0, 0);
            for (int d = 0; d < 1; ++d)
            {
                internalEAreaBegin[d] = std::max(updateEAreaBegin[d], pml->leftDims[d]);
                internalEAreaEnd[d] = std::min(updateEAreaEnd[d],
                    grid->numCells[d] - pml->rightDims[d]);
            }

            const FP coeffCurrent = -(FP)4 * constants::pi * dt;
            const FP cdt = constants::c * dt;
            const FP coeffXY = cdt / (grid->steps.x * anisotropyCoeff.y);
            const FP coeffXZ = cdt / (grid->steps.x * anisotropyCoeff.z);

            // In internal area use:
            // e.x(i, j, k) += dt * -4pi * j.x(i, j, k) + c * dt * ((b.z(i, j+1, k) -
            //     b.z(i, j, k)) / eps_y * dy - (b.y(i, j, k+1) - b.y(i, j, k)) / eps_z * dz),
            // e.y(i, j, k) += dt * -4pi * j.y(i, j, k) + c * dt * ((b.x(i, j, k+1) -
            //     b.x(i, j, k)) / eps_z * dz - (b.z(i+1, j, k) - b.z(i, j, k)) / eps_x * dx),
            // e.z(i, j, k) += dt * -4pi * j.z(i, j, k) + c * dt * ((b.y(i+1, j, k) -
            //     b.y(i, j, k)) / eps_x * dx - (b.x(i, j+1, k) - b.x(i, j, k)) / eps_y * dy),
            try {

                sycl::id<1> start(internalEAreaBegin.x);
                sycl::range<1> end(internalEAreaEnd.x);

                auto kernel = [&](sycl::handler& h) {

                    auto By = bufBy.get_access<sycl::access::mode::read>(h);
                    auto Bz = bufBz.get_access<sycl::access::mode::read>(h);
                    auto Jx = bufJx.get_access<sycl::access::mode::read>(h);
                    auto Jy = bufJy.get_access<sycl::access::mode::read>(h);
                    auto Jz = bufJz.get_access<sycl::access::mode::read>(h);

                    auto Ex = bufEx.get_access<sycl::access::mode::read_write>(h);
                    auto Ey = bufEy.get_access<sycl::access::mode::read_write>(h);
                    auto Ez = bufEz.get_access<sycl::access::mode::read_write>(h);

                    h.parallel_for(end, start, [=](sycl::item<1> ind) {
                        size_t i = ind[0], j = 0, k = 0;
                        Ex[i][j][k] += coeffCurrent * Jx[i][j][k];
                        Ey[i][j][k] += coeffCurrent * Jy[i][j][k] -
                            coeffXY * (Bz[i + 1][j][k] - Bz[i][j][k]);
                        Ez[i][j][k] += coeffCurrent * Jz[i][j][k] +
                            coeffXZ * (By[i + 1][j][k] - By[i][j][k]);
                    });
                };

                if (device == Devices::CPU)
                    node.cpu_device.submit(kernel).wait_and_throw();
                if (device == Devices::GPU)
                    node.gpu_device.submit(kernel).wait_and_throw();

                {//update data
                    bufEx.get_access<sycl::access::mode::read_write>();
                    bufEy.get_access<sycl::access::mode::read_write>();
                    bufEz.get_access<sycl::access::mode::read_write>();
                }
            }
            catch (sycl::exception& e)
            {
                std::cout << e.what() << std::endl;
            }
            catch (std::exception& e)
            {
                std::cout << e.what() << std::endl;
            }
        }
    }
}
