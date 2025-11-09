#include "../include/BinaryAlloyCluster.h"

BinaryAlloyCluster::BinaryAlloyCluster(int n, const std::string& elemA, const std::string& elemB)
    : numAtoms(n), elementA(elemA), elementB(elemB), energy(0.0) {
    coordinates.resize(3 * n, 0.0);
    atomTypes.resize(n, 1);
}

BinaryAlloyCluster::BinaryAlloyCluster(int nA, int nB, const std::string& elemA, const std::string& elemB)
    : numAtoms(nA + nB), elementA(elemA), elementB(elemB), energy(0.0) {
    coordinates.resize(3 * numAtoms, 0.0);
    atomTypes.resize(numAtoms);
    for (int i = 0; i < nA; ++i) atomTypes[i] = 0;
    for (int i = nA; i < numAtoms; ++i) atomTypes[i] = 1;
}

void BinaryAlloyCluster::setElements(const std::string& elemA, const std::string& elemB) {
    elementA = elemA;
    elementB = elemB;
}

std::string BinaryAlloyCluster::getElementSymbol(int atomIndex) const {
    return (atomTypes[atomIndex] == 0) ? elementA : elementB;
}

std::string BinaryAlloyCluster::getCompositionString() const {
    int nA = getNumElementA();
    int nB = getNumElementB();
    return elementA + std::to_string(nA) + elementB + std::to_string(nB);
}

std::array<double, 3> BinaryAlloyCluster::getAtomPosition(int i) const {
    return { coordinates[i], coordinates[numAtoms + i], coordinates[2 * numAtoms + i] };
}

void BinaryAlloyCluster::setAtomPosition(int i, double x, double y, double z) {
    coordinates[i] = x;
    coordinates[numAtoms + i] = y;
    coordinates[2 * numAtoms + i] = z;
}

int BinaryAlloyCluster::getNumElementA() const {
    return std::count(atomTypes.begin(), atomTypes.end(), 0);
}

int BinaryAlloyCluster::getNumElementB() const {
    return std::count(atomTypes.begin(), atomTypes.end(), 1);
}

void BinaryAlloyCluster::randomInitialize(double boxSize) {
    double r0 = boxSize * std::pow(static_cast<double>(numAtoms), 1.0 / 3.0);
    
    for (size_t i = 0; i < coordinates.size(); ++i) {
            coordinates[i] = (RandomGenerator::uniform() - 0.5) * r0;
    }
}

int BinaryAlloyCluster::countStructureFiles() {
    std::string directory = "data/initial_structures";
    int count = 0;

    if (!std::filesystem::exists(directory)) {
        std::cerr << "Directory not found" << std::endl;
        return 0;
    }

    for (const auto& entry : std::filesystem::directory_iterator(directory)) {
        if (entry.is_regular_file() && entry.path().extension() == ".txt") {
            count++;
        }
    }

    return count;
}

bool BinaryAlloyCluster::loadStructureInitialize(int x,int expectedAtoms)  {
    std::string directory = "data/initial_structures";
    std::string filename = "tmp" + std::to_string(x) + ".txt";
    std::string filepath = directory + "/" + filename;

    if (!std::filesystem::exists(filepath)) {
        std::cerr << "File not found: " << filename << std::endl;
        return false;
    }

    std::ifstream file(filepath);
    if (!file.is_open()) {
        std::cerr << "Cannot open: " << filename << std::endl;
        return false;
    }

    std::vector<std::array<double, 3>> positions;
    std::string line;

    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;

        std::istringstream iss(line);
        std::string symbol;
        double x_coord, y_coord, z_coord;

        if (iss >> symbol >> x_coord >> y_coord >> z_coord) {
            positions.push_back({ x_coord, y_coord, z_coord });
        }
    }

    if (positions.size() != expectedAtoms) {
        std::cerr << "Error: " << filename << " has " << positions.size()
            << " atoms, expected " << expectedAtoms << std::endl;
        return false;
    }

    //std::cout << filename << ": " << positions.size() << " atoms (OK)" << std::endl;
    for (size_t i = 0; i < coordinates.size()/3; ++i) {
        coordinates[i] = positions[i][0];                // x
        coordinates[expectedAtoms + i] = positions[i][1];     // y
        coordinates[2 * expectedAtoms + i] = positions[i][2]; // z
    }
    return true;
}

void BinaryAlloyCluster::centerAtOrigin() {
    double cx = 0, cy = 0, cz = 0;
    for (int i = 0; i < numAtoms; ++i) {
        cx += coordinates[i];
        cy += coordinates[numAtoms + i];
        cz += coordinates[2 * numAtoms + i];
    }
    cx /= numAtoms;
    cy /= numAtoms;
    cz /= numAtoms;
    
    for (int i = 0; i < numAtoms; ++i) {
        coordinates[i] -= cx;
        coordinates[numAtoms + i] -= cy;
        coordinates[2 * numAtoms + i] -= cz;
    }
}

double BinaryAlloyCluster::getRadius() const {
    double maxDist = 0.0;
    for (int i = 0; i < numAtoms; ++i) {
        double r2 = 0.0;
        r2 += coordinates[i] * coordinates[i];
        r2 += coordinates[numAtoms + i] * coordinates[numAtoms + i];
        r2 += coordinates[2 * numAtoms + i] * coordinates[2 * numAtoms + i];
        maxDist = (std::max)(maxDist, r2);
    }
    return std::sqrt(maxDist);
}
