#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <omp.h>

void precomputeTrigonometricValues(int latitudeDivisions, int longitudeDivisions,
                                   std::vector<float>& sinTheta, std::vector<float>& cosTheta,
                                   std::vector<float>& sinPhi, std::vector<float>& cosPhi) {
    for (int i = 0; i <= latitudeDivisions; ++i) {
        float theta = i * M_PI / latitudeDivisions;
        sinTheta[i] = sin(theta);
        cosTheta[i] = cos(theta);
    }
    for (int j = 0; j <= longitudeDivisions; ++j) {
        float phi = j * 2 * M_PI / longitudeDivisions;
        sinPhi[j] = sin(phi);
        cosPhi[j] = cos(phi);
    }
}

void generateSphereFrame(std::vector<float>& frameData, float radius, int latitudeDivisions, int longitudeDivisions,
                         const std::vector<float>& sinTheta, const std::vector<float>& cosTheta,
                         const std::vector<float>& sinPhi, const std::vector<float>& cosPhi,
                         float deformationFactor, float heightFactor, float time, int n, int m) {
    #pragma omp parallel for collapse(2)
    for (int i = 0; i <= latitudeDivisions; ++i) {
        for (int j = 0; j <= longitudeDivisions; ++j) {
            int index = (i * (longitudeDivisions + 1) + j) * 3;

            float deformation = 1.0f - deformationFactor * cosTheta[i] * cosTheta[i];
            float x = radius * sinTheta[i] * cosPhi[j];
            float y = radius * sinTheta[i] * sinPhi[j];
            float z = radius * cosTheta[i] * deformation * heightFactor;

            float oscillation = 0.1f * sin(n * acos(cosTheta[i])) * cos(m * atan2(sinPhi[j], cosPhi[j]) + time);
            z += oscillation;

            frameData[index] = x;
            frameData[index + 1] = y;
            frameData[index + 2] = z;
        }
    }
}

void saveFrameToFile(std::ofstream& outFile, const std::vector<float>& frameData, int latitudeDivisions, int longitudeDivisions) {
    for (int i = 0; i <= latitudeDivisions; ++i) {
        for (int j = 0; j <= longitudeDivisions; ++j) {
            int index = (i * (longitudeDivisions + 1) + j) * 3;
            outFile << frameData[index] << " " << frameData[index + 1] << " " << frameData[index + 2] << "\n";
        }
        outFile << "\n"; // Separate rows for easier parsing
    }
}

int main() {
    float radius = 1.0f;               // Radius of the sphere
    int latitudeDivisions = 200;      // Higher resolution
    int longitudeDivisions = 200;     // Higher resolution
    float deformationFactor = 0.8f;   // Adjustable concave deformation factor
    float heightFactor = 0.8f;        // Adjustable height factor
    int totalFrames = 1000;            // Number of frames for animation
    float timeStep = 0.1f;            // Time increment per frame
    int n = 5;                        // Energy level n (number of oscillations in theta direction)
    int m = 5;                        // Energy level m (number of oscillations in phi direction)
    std::string outputFile = "red_blood_cell_dynamics_data.txt";

    // Allocate precomputed trigonometric values
    std::vector<float> sinTheta(latitudeDivisions + 1), cosTheta(latitudeDivisions + 1);
    std::vector<float> sinPhi(longitudeDivisions + 1), cosPhi(longitudeDivisions + 1);

    // Precompute trigonometric values
    precomputeTrigonometricValues(latitudeDivisions, longitudeDivisions, sinTheta, cosTheta, sinPhi, cosPhi);

    // Allocate memory for a single frame
    std::vector<float> frameData((latitudeDivisions + 1) * (longitudeDivisions + 1) * 3);

    // Open the output file for all frames
    std::ofstream outFile(outputFile);
    if (!outFile) {
        std::cerr << "Error opening file for writing: " << outputFile << std::endl;
        return 1;
    }

    // Write the frames
    for (int frame = 0; frame < totalFrames; ++frame) {
        float time = frame * timeStep;

        // Generate the frame
        generateSphereFrame(frameData, radius, latitudeDivisions, longitudeDivisions,
                            sinTheta, cosTheta, sinPhi, cosPhi,
                            deformationFactor, heightFactor, time, n, m);

        // Write to the file
        outFile << "# Frame " << frame << " (n=" << n << ", m=" << m << ")\n";
        saveFrameToFile(outFile, frameData, latitudeDivisions, longitudeDivisions);
    }

    outFile.close();
    std::cout << "All animation frames with energy levels saved to " << outputFile << std::endl;

    return 0;
}
