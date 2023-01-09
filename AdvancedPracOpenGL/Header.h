// Functionality to compare data points and return the Euclidean distance between them

#pragma once
#include <vector>

// Return the Euclidean distance between two vectors
float euclideanDistance(std::vector<float>* vec1, std::vector<float>* vec2)
{
	float sum = 0;
	std::vector<float>::iterator it2 = vec2->begin();

	// Sum the squares of the differences of the matching attributes between vectors
	for (auto& attr : *vec1)
	{
		sum += pow(attr - *it2, 2);
	}

	return sqrt(sum);
}

