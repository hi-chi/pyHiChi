#pragma once
#include "pyFieldInterface.h"
#include "ScalarField.h"

#include "pybind11/pybind11.h"

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

        virtual FP3 getE(const FP3& coords) const = 0;
        virtual FP3 getB(const FP3& coords) const = 0;
        virtual FP3 getJ(const FP3& coords) const = 0;

        FP3 getE(FP x, FP y, FP z) const { return getE(FP3(x, y, z)); }
        FP3 getB(FP x, FP y, FP z) const { return getB(FP3(x, y, z)); }
        FP3 getJ(FP x, FP y, FP z) const { return getJ(FP3(x, y, z)); }

        void getFields(const FP3& coords, FP3& e, FP3& b) const {
            e = getE(coords);
            b = getB(coords);
        }

        virtual void updateFields() = 0;
        virtual void advance(FP dt) = 0;

        virtual std::shared_ptr<pyFieldBase> applyMapping(
            const std::shared_ptr<pyFieldBase>& self,
            const std::shared_ptr<Mapping>& mapping) const = 0;
    };


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

        std::shared_ptr<pyScalarField> getEx() {
            return std::make_shared<pyScalarField>(&(this->getField()->getGrid()->Ex));
        }
        std::shared_ptr<pyScalarField> getEy() {
            return std::make_shared<pyScalarField>(&(this->getField()->getGrid()->Ey));
        }
        std::shared_ptr<pyScalarField> getEz() {
            return std::make_shared<pyScalarField>(&(this->getField()->getGrid()->Ez));
        }

        std::shared_ptr<pyScalarField> getBx() {
            return std::make_shared<pyScalarField>(&(this->getField()->getGrid()->Bx));
        }
        std::shared_ptr<pyScalarField> getBy() {
            return std::make_shared<pyScalarField>(&(this->getField()->getGrid()->By));
        }
        std::shared_ptr<pyScalarField> getBz() {
            return std::make_shared<pyScalarField>(&(this->getField()->getGrid()->Bz));
        }

        std::shared_ptr<pyScalarField> getJx() {
            return std::make_shared<pyScalarField>(&(this->getField()->getGrid()->Jx));
        }
        std::shared_ptr<pyScalarField> getJy() {
            return std::make_shared<pyScalarField>(&(this->getField()->getGrid()->Jy));
        }
        std::shared_ptr<pyScalarField> getJz() {
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
        {
            pyFieldPointer =
                std::dynamic_pointer_cast<pyField<TGrid, TFieldSolver>>(pyWrappedField).get();
            pyMappedFieldPointer =
                std::dynamic_pointer_cast<pyMappedField<TGrid, TFieldSolver>>(pyWrappedField).get();
            if (!pyFieldPointer && !pyMappedFieldPointer)
                std::cout << "ERROR in pyMappedField::getField(): wrong cast" << std::endl;
        }

        inline Field<TGrid, TFieldSolver>* getField() const {
            if (pyFieldPointer) return pyFieldPointer->getField();
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
            bool status = true;
            FP time = getField()->getFieldSolver()->globalTime +
                getField()->getFieldSolver()->timeShiftE;
            FP3 inverseCoords = getInverseCoords(coords, time, &status);
            if (!status) return FP3(0, 0, 0);
            return BaseInterface::getE(inverseCoords);
        }

        FP3 getB(const FP3& coords) const override {
            bool status = true;
            FP time = getField()->getFieldSolver()->globalTime +
                getField()->getFieldSolver()->timeShiftB;
            FP3 inverseCoords = getInverseCoords(coords, time, &status);
            if (!status) return FP3(0, 0, 0);
            return BaseInterface::getB(inverseCoords);
        }

        FP3 getJ(const FP3& coords) const override {
            bool status = true;
            FP time = getField()->getFieldSolver()->globalTime +
                getField()->getFieldSolver()->timeShiftJ;
            FP3 inverseCoords = getInverseCoords(coords, time, &status);
            if (!status) return FP3(0, 0, 0);
            return BaseInterface::getJ(inverseCoords);
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

        // these pointers are used not to do dynamic cast each time
        // always one of them is NULL
        pyField<TGrid, TFieldSolver>* pyFieldPointer;
        pyMappedField<TGrid, TFieldSolver>* pyMappedFieldPointer;

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
