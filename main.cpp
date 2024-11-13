#include <iostream>
#include <cmath>

const double epsilon_0 = 8.854187817e-12;
const double mu_0 = 4 * M_PI * 1e-7;

class Field {
protected:
    double *value;

public:
    Field() : value(new double[3]{0.0, 0.0, 0.0}) {}

    Field(double x, double y, double z) : value(new double[3]{x, y, z}) {}

    virtual ~Field() {
        delete[] value;
    }

    Field(const Field &other) {
        value = new double[3];
        for (int i = 0; i < 3; ++i) {
            value[i] = other.value[i];
        }
    }

    virtual void printMagnitude() const {
        std::cout << "Field components:(" << value[0] << ", " << value[1] << ", " << value[2] << ")\n";
    }
    friend std::ostream &operator<<(std::ostream &os, const Field &field) {
        os << "(" << field.value[0] << ", " << field.value[1] << ", " << field.value[2] << ")";
        return os;
    }
};
class ElectricField : public Field {
private:
    double calculatedField;

public:
    ElectricField() : Field(), calculatedField(0.0) {}

    ElectricField(double x, double y, double z) : Field(x, y, z), calculatedField(0.0) {}

    void calculateElectricField(double Q, double r) {
        calculatedField = Q / (4 * M_PI * epsilon_0 * r * r);
    }
    ElectricField operator+(const ElectricField &other) const {
        return ElectricField(value[0] + other.value[0], value[1] + other.value[1], value[2] + other.value[2]);
    }

    void printCalculatedField() const {
        std::cout << "Calculated Electric Field : " << calculatedField << " N/C\n";
    }
};
class MagneticField : public Field {
private:
double calculatedField;

public:
    MagneticField() : Field(), calculatedField(0.0) {}

    MagneticField(double x, double y, double z) : Field(x, y, z), calculatedField(0.0) {}

    void calculateMagneticField(double I, double r) {
        calculatedField = (mu_0 * I) / (2 * M_PI * r);
    }

    MagneticField operator+(const MagneticField &other) const {
        return MagneticField(value[0] + other.value[0], value[1] + other.value[1], value[2] + other.value[2]);
    }

    void printCalculatedField() const {
        std::cout << "Calculated Magnetic Field : " << calculatedField << " T\n";
    }
};

int main() {
    ElectricField eField1(1.0, 2.0, 3.0);
    MagneticField mField1(0.5, 1.0, 1.5);

    std::cout << "Initial Electric Field Components: " << eField1 << std::endl;
    std::cout << "Initial Magnetic Field Components: " << mField1 << std::endl;

    double charge = 1e-6;
    double distance_electric = 0.1;
    eField1.calculateElectricField(charge, distance_electric);
    std::cout << "Electric Field Calculation for Q =  " << charge << " C, r = " << distance_electric << " m\n";
    eField1.printCalculatedField();

    double current = 2.0;
    double distance_magnetic = 0.05;
    mField1.calculateMagneticField(current, distance_magnetic);
    std::cout << "Magnetic Field Calculation for I =  " << current << " A, r = " << distance_magnetic << " m\n";
    mField1.printCalculatedField();

    ElectricField eField2(4.0, 5.0, 6.0);
    std::cout << "Second Electric Field Components: " << eField2 << std::endl;
    ElectricField eField3 = eField1 + eField2;
    std::cout << "Combined Electric Field Components (: " << eField3 << std::endl;
    MagneticField mField2(1.0, 1.5, 2.0);
    std::cout << "Second Magnetic Field Components: " << mField2 << std::endl;
    MagneticField mField3 = mField1 + mField2;
    std::cout << "Combined Magnetic Field Components : " << mField3 << std::endl;

    return 0;
}

