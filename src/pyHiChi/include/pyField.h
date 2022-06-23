#pragma once
#include "pyFieldInterface.h"
#include "ScalarField.h"

#include "pybind11/pybind11.h"
#include "pybind11/numpy.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace pfc
{
    // wrapper over ScalarField class object
    class pyScalarField {
    public:

        pyScalarField(ScalarField<FP>* scalarField) :scalarField(scalarField) {}

        FP* getData() const {
            return scalarField->getData();
        }

        Int3 getSize() const {
            return scalarField->getSize();
        }

        // read-only accessors
        FP get(const Int3& index) const {
            return (*scalarField)(index);
        }

        FP get(int i, int j, int k) const {
            return (*scalarField)(i, j, k);
        }

    private:

        ScalarField<FP>* scalarField;
    };


    // Base class for all pyFields
    class pyFieldBase {
    public:

        virtual FP getEx(const FP3& coords) const = 0;
        virtual FP getEy(const FP3& coords) const = 0;
        virtual FP getEz(const FP3& coords) const = 0;

        virtual FP getBx(const FP3& coords) const = 0;
        virtual FP getBy(const FP3& coords) const = 0;
        virtual FP getBz(const FP3& coords) const = 0;

        virtual FP getJx(const FP3& coords) const = 0;
        virtual FP getJy(const FP3& coords) const = 0;
        virtual FP getJz(const FP3& coords) const = 0;

        virtual FP3 getE(const FP3& coords) const = 0;
        virtual FP3 getB(const FP3& coords) const = 0;
        virtual FP3 getJ(const FP3& coords) const = 0;

        void getFields(const FP3& coords, FP3& e, FP3& b) const {
            e = getE(coords);
            b = getB(coords);
        }

        virtual void updateFields() = 0;
        virtual void advance(FP dt) = 0;

        virtual std::shared_ptr<pyFieldBase> applyMapping(
            const std::shared_ptr<pyFieldBase>& self,
            const std::shared_ptr<Mapping>& mapping) const = 0;

        // functions to get cross sections
        py::array_t<FP> getSlice1d(
            CoordinateEnum crossAxis1, FP pos1,
            CoordinateEnum crossAxis2, FP pos2,
            CoordinateEnum axis, FP minCoord, FP maxCoord, size_t size,
            FP(pyFieldBase::*getFieldValue)(const FP3&) const)
        {
            py::array_t<FP> res({ static_cast<pybind11::size_t>(size) });
            auto accRes = res.mutable_unchecked<1>();
            FP step = (maxCoord - minCoord) / size;
            OMP_FOR()
            for (py::ssize_t i = 0; i < size; i++) {
                FP3 coords;
                coords[(int)crossAxis1] = pos1;
                coords[(int)crossAxis1] = pos2;
                coords[(int)axis] = minCoord + step * i;
                accRes(i) = (this->*getFieldValue)(coords);
            }
            return res;
        }

        py::array_t<FP> getSlice2d(
            CoordinateEnum crossAxis, FP pos,
            CoordinateEnum axis1, FP minCoord1, FP maxCoord1, size_t size1,
            CoordinateEnum axis2, FP minCoord2, FP maxCoord2, size_t size2,
            FP(pyFieldBase::*getFieldValue)(const FP3&) const)
        {
            py::array_t<FP> res({ static_cast<pybind11::size_t>(size1),
                                  static_cast<pybind11::size_t>(size2) });
            auto accRes = res.mutable_unchecked<2>();
            FP step1 = (maxCoord1 - minCoord1) / size1,
                step2 = (maxCoord2 - minCoord2) / size2;
            OMP_FOR()
            for (py::ssize_t i = 0; i < size1; i++)
                for (py::ssize_t j = 0; j < size2; j++) {
                    FP3 coords;
                    coords[(int)crossAxis] = pos;
                    coords[(int)axis1] = minCoord1 + step1 * i;
                    coords[(int)axis2] = minCoord2 + step2 * j;
                    accRes(i, j) = (this->*getFieldValue)(coords);
                }
            return res;
        }

        py::array_t<FP> getSlice3d(
            CoordinateEnum axis1, FP minCoord1, FP maxCoord1, size_t size1,
            CoordinateEnum axis2, FP minCoord2, FP maxCoord2, size_t size2,
            CoordinateEnum axis3, FP minCoord3, FP maxCoord3, size_t size3,
            FP(pyFieldBase::*getFieldValue)(const FP3&) const)
        {
            py::array_t<FP> res({ static_cast<pybind11::size_t>(size1),
                                  static_cast<pybind11::size_t>(size2),
                                  static_cast<pybind11::size_t>(size3) });
            auto accRes = res.mutable_unchecked<3>();
            FP step1 = (maxCoord1 - minCoord1) / size1,
                step2 = (maxCoord2 - minCoord2) / size2,
                step3 = (maxCoord3 - minCoord3) / size3;
            OMP_FOR()
            for (py::ssize_t i = 0; i < size1; i++)
                for (py::ssize_t j = 0; j < size2; j++)
                    for (py::ssize_t k = 0; k < size3; k++) {
                        FP3 coords;
                        coords[(int)axis1] = minCoord1 + step1 * i;
                        coords[(int)axis2] = minCoord2 + step2 * j;
                        coords[(int)axis3] = minCoord3 + step3 * k;
                        accRes(i, j, k) = (this->*getFieldValue)(coords);
                    }
            return res;
        }
    };

    template<class TGrid, class TFieldSolver>
    class pyMappedField;

    // Simple pyField class
    template<class TGrid, class TFieldSolver>
    class pyField : public pyFieldInterface<TGrid, TFieldSolver, pyField<TGrid, TFieldSolver>>,
        public pyFieldBase
    {
        using BaseInterface = pyFieldInterface<TGrid, TFieldSolver, pyField<TGrid, TFieldSolver>>;
    
    public:
    
        pyField(const Int3 & numInternalCells,
            const FP3 & minCoords, const FP3 & steps, FP dt) :
            field(new Field<TGrid, TFieldSolver>(numInternalCells,
                minCoords, steps, dt))
        {}
    
        pyField(FP dt) :
            field(new Field<TGrid, TFieldSolver>(dt))
        {}
    
        inline Field<TGrid, TFieldSolver>* getField() const {
            return field.get();
        }
    
        inline FP3 convertCoords(const FP3& coords, FP timeShift = 0.0) const {
            return coords;
        }
    
        std::shared_ptr<pyFieldBase> applyMapping(
            const std::shared_ptr<pyFieldBase>& self,
            const std::shared_ptr<Mapping>& mapping) const override {
            return std::static_pointer_cast<pyFieldBase>(
                std::make_shared<pyMappedField<TGrid, TFieldSolver>>(
                    std::static_pointer_cast<pyFieldBase>(self), mapping
                    )
                );
        }
    
        FP3 getE(const FP3& coords) const override {
            return BaseInterface::getE(coords);
        }
        FP3 getB(const FP3& coords) const override {
            return BaseInterface::getB(coords);
        }
        FP3 getJ(const FP3& coords) const override {
            return BaseInterface::getJ(coords);
        }

        FP getEx(const FP3& coords) const override {
            return BaseInterface::getEx(coords);
        }
        FP getEy(const FP3& coords) const override {
            return BaseInterface::getEy(coords);
        }
        FP getEz(const FP3& coords) const override {
            return BaseInterface::getEz(coords);
        }

        FP getBx(const FP3& coords) const override {
            return BaseInterface::getBx(coords);
        }
        FP getBy(const FP3& coords) const override {
            return BaseInterface::getBy(coords);
        }
        FP getBz(const FP3& coords) const override {
            return BaseInterface::getBz(coords);
        }

        FP getJx(const FP3& coords) const override {
            return BaseInterface::getJx(coords);
        }
        FP getJy(const FP3& coords) const override {
            return BaseInterface::getJy(coords);
        }
        FP getJz(const FP3& coords) const override {
            return BaseInterface::getJz(coords);
        }
    
        void updateFields() override {
            return BaseInterface::updateFields();
        }
    
        void advance(FP dt) override {
            return BaseInterface::advance(dt);
        }

        std::shared_ptr<pyField<TGrid, TFieldSolver>> zoom(const FP3& minCoord,
            const FP3& zoomedGridSize, const FP3& zoomedGridStep) const {
            std::shared_ptr<pyField<TGrid, TFieldSolver>> zoomedField;
            zoomedField.reset(new pyField<TGrid, TFieldSolver>(
                (Int3)zoomedGridSize, minCoord, zoomedGridStep,
                this->getField()->getFieldSolver()->dt)
            );
            this->getField()->getGrid()->copyValues(zoomedField->getField()->getGrid());
            return zoomedField;
        }

        std::shared_ptr<pyScalarField> getExArray() {
            return std::make_shared<pyScalarField>(&(this->getField()->getGrid()->Ex));
        }
        std::shared_ptr<pyScalarField> getEyArray() {
            return std::make_shared<pyScalarField>(&(this->getField()->getGrid()->Ey));
        }
        std::shared_ptr<pyScalarField> getEzArray() {
            return std::make_shared<pyScalarField>(&(this->getField()->getGrid()->Ez));
        }

        std::shared_ptr<pyScalarField> getBxArray() {
            return std::make_shared<pyScalarField>(&(this->getField()->getGrid()->Bx));
        }
        std::shared_ptr<pyScalarField> getByArray() {
            return std::make_shared<pyScalarField>(&(this->getField()->getGrid()->By));
        }
        std::shared_ptr<pyScalarField> getBzArray() {
            return std::make_shared<pyScalarField>(&(this->getField()->getGrid()->Bz));
        }

        std::shared_ptr<pyScalarField> getJxArray() {
            return std::make_shared<pyScalarField>(&(this->getField()->getGrid()->Jx));
        }
        std::shared_ptr<pyScalarField> getJyArray() {
            return std::make_shared<pyScalarField>(&(this->getField()->getGrid()->Jy));
        }
        std::shared_ptr<pyScalarField> getJzArray() {
            return std::make_shared<pyScalarField>(&(this->getField()->getGrid()->Jz));
        }
    
    private:
    
        std::unique_ptr<Field<TGrid, TFieldSolver>> field;
    
    };

    typedef pyField<YeeGrid, FDTD> pyYeeField;
    typedef pyField<PSTDGrid, PSTD> pyPSTDField;
    typedef pyField<PSATDGrid, PSATD> pyPSATDField;
    typedef pyField<PSATDGrid, PSATDPoisson> pyPSATDPoissonField;
    typedef pyField<PSATDTimeStraggeredGrid, PSATDTimeStraggered> pyPSATDTimeStraggeredField;
    typedef pyField<PSATDTimeStraggeredGrid, PSATDTimeStraggeredPoisson> pyPSATDTimeStraggeredPoissonField;

    typedef pyField<NoGrid, NoFieldSolver> pyAnalyticalField;

    template<>
    pyAnalyticalField::pyField(FP dt) :
        field(new Field<NoGrid, NoFieldSolver>(dt)) {}


    // pyField class that supports mappings
    // wrapper over pyField class object or pyMappedField class object
    template<class TGrid, class TFieldSolver>
    class pyMappedField : public pyFieldInterface<TGrid, TFieldSolver,
        pyMappedField<TGrid, TFieldSolver>>, public pyFieldBase
    {
        using BaseInterface = pyFieldInterface<TGrid, TFieldSolver, pyMappedField<TGrid, TFieldSolver>>;

    public:

        pyMappedField(const std::shared_ptr<pyFieldBase>& other,
            const std::shared_ptr<Mapping>& mapping) :
            pyWrappedField(other), mapping(mapping)
        {}

        inline Field<TGrid, TFieldSolver>* getField() const {
            std::shared_ptr<pyField<TGrid, TFieldSolver>> pyFieldPointer =
                std::dynamic_pointer_cast<pyField<TGrid, TFieldSolver>>(pyWrappedField);
            if (pyFieldPointer) return pyFieldPointer->getField();

            std::shared_ptr<pyMappedField<TGrid, TFieldSolver>> pyMappedFieldPointer =
                std::dynamic_pointer_cast<pyMappedField<TGrid, TFieldSolver>>(pyWrappedField);
            if (pyMappedFieldPointer) return pyMappedFieldPointer->getField();

            return 0;
        }

        inline FP3 convertCoords(const FP3& coords, FP timeShift = 0.0) const {
            bool status = true;
            return getDirectCoords(coords,
                getField()->getFieldSolver()->globalTime + timeShift, &status);
        }

        std::shared_ptr<pyFieldBase> applyMapping(
            const std::shared_ptr<pyFieldBase>& self,
            const std::shared_ptr<Mapping>& mapping) const override {
            return std::static_pointer_cast<pyFieldBase>(
                std::make_shared<pyMappedField<TGrid, TFieldSolver>>(
                    std::static_pointer_cast<pyFieldBase>(self), mapping
                    )
                );
        }

        FP3 getE(const FP3& coords) const override {
            return this->getFieldComp3(coords, getField()->getFieldSolver()->timeShiftE,
                &BaseInterface::getE);
        }
        FP3 getB(const FP3& coords) const override {
            return this->getFieldComp3(coords, getField()->getFieldSolver()->timeShiftB,
                &BaseInterface::getB);
        }
        FP3 getJ(const FP3& coords) const override {
            return this->getFieldComp3(coords, getField()->getFieldSolver()->timeShiftJ,
                &BaseInterface::getJ);
        }

        FP getEx(const FP3& coords) const override {
            return this->getFieldComp(coords, getField()->getFieldSolver()->timeShiftE,
                &BaseInterface::getEx);
        }
        FP getEy(const FP3& coords) const override {
            return this->getFieldComp(coords, getField()->getFieldSolver()->timeShiftE,
                &BaseInterface::getEy);
        }
        FP getEz(const FP3& coords) const override {
            return this->getFieldComp(coords, getField()->getFieldSolver()->timeShiftE,
                &BaseInterface::getEz);
        }

        FP getBx(const FP3& coords) const override {
            return this->getFieldComp(coords, getField()->getFieldSolver()->timeShiftB,
                &BaseInterface::getBx);
        }
        FP getBy(const FP3& coords) const override {
            return this->getFieldComp(coords, getField()->getFieldSolver()->timeShiftB,
                &BaseInterface::getBy);
        }
        FP getBz(const FP3& coords) const override {
            return this->getFieldComp(coords, getField()->getFieldSolver()->timeShiftB,
                &BaseInterface::getBz);
        }

        FP getJx(const FP3& coords) const override {
            return this->getFieldComp(coords, getField()->getFieldSolver()->timeShiftJ,
                &BaseInterface::getJx);
        }
        FP getJy(const FP3& coords) const override {
            return this->getFieldComp(coords, getField()->getFieldSolver()->timeShiftJ,
                &BaseInterface::getJy);
        }
        FP getJz(const FP3& coords) const override {
            return this->getFieldComp(coords, getField()->getFieldSolver()->timeShiftJ,
                &BaseInterface::getJz);
        }

        void updateFields() override {
            return BaseInterface::updateFields();
        }

        void advance(FP dt) override {
            return BaseInterface::advance(dt);
        }

    protected:

        inline FP3 getDirectCoords(const FP3& coords, FP time, bool* status) const {
            FP3 coords_ = coords;
            *status = true;
            std::shared_ptr<pyMappedField<TGrid, TFieldSolver>> pyMappedFieldPointer =
                std::dynamic_pointer_cast<pyMappedField<TGrid, TFieldSolver>>(pyWrappedField);
            if (pyMappedFieldPointer)
                coords_ = pyMappedFieldPointer->getDirectCoords(coords_, time, status);
            bool status2 = true;
            if (mapping) coords_ = mapping->getDirectCoords(coords_, time, &status2);
            *status = *status && status2;
            return coords_;
        }

        inline FP3 getInverseCoords(const FP3& coords, FP time, bool* status) const {
            FP3 coords_ = coords;
            *status = true;
            std::shared_ptr<pyMappedField<TGrid, TFieldSolver>> pyMappedFieldPointer =
                std::dynamic_pointer_cast<pyMappedField<TGrid, TFieldSolver>>(pyWrappedField);
            if (pyMappedFieldPointer)
                coords_ = pyMappedFieldPointer->getInverseCoords(coords_, time, status);
            bool status2 = true;
            if (mapping) coords_ = mapping->getInverseCoords(coords_, time, &status2);
            *status = *status && status2;
            return coords_;
        }

    private:

        std::shared_ptr<pyFieldBase> pyWrappedField;
        std::shared_ptr<Mapping> mapping;

        FP getFieldComp(const FP3& coords, FP timeShift,
            FP(BaseInterface::* getFieldValue)(const FP3&) const) const
        {
            bool status = true;
            FP time = getField()->getFieldSolver()->globalTime + timeShift;
            FP3 inverseCoords = getInverseCoords(coords, time, &status);
            if (!status) return 0.0;
            return (this->*getFieldValue)(inverseCoords);
        }

        FP3 getFieldComp3(const FP3& coords, FP timeShift,
            FP3(BaseInterface::* getFieldValue)(const FP3&) const) const
        {
            bool status = true;
            FP time = getField()->getFieldSolver()->globalTime + timeShift;
            FP3 inverseCoords = getInverseCoords(coords, time, &status);
            if (!status) return FP3(0.0, 0.0, 0.0);
            return (this->*getFieldValue)(inverseCoords);
        }

    };

    typedef pyMappedField<YeeGrid, FDTD> pyMappedYeeField;
    typedef pyMappedField<PSTDGrid, PSTD> pyMappedPSTDField;
    typedef pyMappedField<PSATDGrid, PSATD> pyMappedPSATDField;
    typedef pyMappedField<PSATDGrid, PSATDPoisson> pyMappedPSATDPoissonField;
    typedef pyMappedField<PSATDTimeStraggeredGrid, PSATDTimeStraggered> pyMappedPSATDTimeStraggeredField;
    typedef pyMappedField<PSATDTimeStraggeredGrid, PSATDTimeStraggeredPoisson> pyMappedPSATDTimeStraggeredPoissonField;

    typedef pyMappedField<NoGrid, NoFieldSolver> pyMappedAnalyticalField;


    // Object returned when summing fields
    class pySumField : public pyFieldBase
    {
    public:

        pySumField(const std::shared_ptr<pyFieldBase>& pyWrappedField1,
            const std::shared_ptr<pyFieldBase>& pyWrappedField2) :
            pyWrappedField1(pyWrappedField1), pyWrappedField2(pyWrappedField2)
        {}

        pySumField(const std::shared_ptr<pySumField>& other,
            const std::shared_ptr<Mapping>& mapping) :
            pyWrappedField1(other->pyWrappedField1->applyMapping(other->pyWrappedField1, mapping)),
            pyWrappedField2(other->pyWrappedField2->applyMapping(other->pyWrappedField2, mapping))
        {}

        std::shared_ptr<pyFieldBase> applyMapping(
            const std::shared_ptr<pyFieldBase>& self,
            const std::shared_ptr<Mapping>& mapping) const override {
            return std::static_pointer_cast<pyFieldBase>(
                std::make_shared<pySumField>(
                    std::static_pointer_cast<pySumField>(self), mapping
                    )
                );
        }

        FP3 getE(const FP3& coords) const override {
            return pyWrappedField1->getE(coords) + pyWrappedField2->getE(coords);
        }
        FP3 getB(const FP3& coords) const override {
            return pyWrappedField1->getB(coords) + pyWrappedField2->getB(coords);
        }
        FP3 getJ(const FP3& coords) const override {
            return pyWrappedField1->getJ(coords) + pyWrappedField2->getJ(coords);
        }

        FP getEx(const FP3& coords) const override {
            return pyWrappedField1->getEx(coords) + pyWrappedField2->getEx(coords);
        }
        FP getEy(const FP3& coords) const override {
            return pyWrappedField1->getEy(coords) + pyWrappedField2->getEy(coords);
        }
        FP getEz(const FP3& coords) const override {
            return pyWrappedField1->getEz(coords) + pyWrappedField2->getEz(coords);
        }

        FP getBx(const FP3& coords) const override {
            return pyWrappedField1->getBx(coords) + pyWrappedField2->getBx(coords);
        }
        FP getBy(const FP3& coords) const override {
            return pyWrappedField1->getBy(coords) + pyWrappedField2->getBy(coords);
        }
        FP getBz(const FP3& coords) const override {
            return pyWrappedField1->getBz(coords) + pyWrappedField2->getBz(coords);
        }

        FP getJx(const FP3& coords) const override {
            return pyWrappedField1->getJx(coords) + pyWrappedField2->getJx(coords);
        }
        FP getJy(const FP3& coords) const override {
            return pyWrappedField1->getJy(coords) + pyWrappedField2->getJy(coords);
        }
        FP getJz(const FP3& coords) const override {
            return pyWrappedField1->getJz(coords) + pyWrappedField2->getJz(coords);
        }

        void updateFields() override {
            pyWrappedField1->updateFields();
            pyWrappedField2->updateFields();
        }

        void advance(FP dt) override {
            pyWrappedField1->advance(dt);
            pyWrappedField2->advance(dt);
        }

    private:

        std::shared_ptr<pyFieldBase> pyWrappedField1;
        std::shared_ptr<pyFieldBase> pyWrappedField2;
    };


    // Object returned when multiplying fields by factor
    class pyMulField : public pyFieldBase {
    public:

        pyMulField(const std::shared_ptr<pyFieldBase>& pyWrappedField, FP factor) :
            pyWrappedField(pyWrappedField),
            factor(factor)
        {}

        pyMulField(const std::shared_ptr<pyMulField>& other,
            const std::shared_ptr<Mapping>& mapping) :
            pyWrappedField(other->pyWrappedField->applyMapping(other->pyWrappedField, mapping))
        {}

        std::shared_ptr<pyFieldBase> applyMapping(
            const std::shared_ptr<pyFieldBase>& self,
            const std::shared_ptr<Mapping>& mapping) const override {
            return std::static_pointer_cast<pyFieldBase>(
                std::make_shared<pyMulField>(
                    std::static_pointer_cast<pyMulField>(self), mapping
                    )
                );
        }

        FP3 getE(const FP3& coords) const override {
            return pyWrappedField->getE(coords) * factor;
        }
        FP3 getB(const FP3& coords) const override {
            return pyWrappedField->getB(coords) * factor;
        }
        FP3 getJ(const FP3& coords) const override {
            return pyWrappedField->getJ(coords) * factor;
        }

        FP getEx(const FP3& coords) const override {
            return pyWrappedField->getEx(coords) * factor;
        }
        FP getEy(const FP3& coords) const override {
            return pyWrappedField->getEy(coords) * factor;
        }
        FP getEz(const FP3& coords) const override {
            return pyWrappedField->getEz(coords) * factor;
        }

        FP getBx(const FP3& coords) const override {
            return pyWrappedField->getBx(coords) * factor;
        }
        FP getBy(const FP3& coords) const override {
            return pyWrappedField->getBy(coords) * factor;
        }
        FP getBz(const FP3& coords) const override {
            return pyWrappedField->getBz(coords) * factor;
        }

        FP getJx(const FP3& coords) const override {
            return pyWrappedField->getJx(coords) * factor;
        }
        FP getJy(const FP3& coords) const override {
            return pyWrappedField->getJy(coords) * factor;
        }
        FP getJz(const FP3& coords) const override {
            return pyWrappedField->getJz(coords) * factor;
        }

        void updateFields() override {
            pyWrappedField->updateFields();
        }

        void advance(FP dt) override {
            pyWrappedField->advance(dt);
        }

    private:

        FP factor = 1.0;
        std::shared_ptr<pyFieldBase> pyWrappedField;
    };

}
