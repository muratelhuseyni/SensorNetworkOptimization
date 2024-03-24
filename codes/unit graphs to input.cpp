//#include "Common.h"
//#include <bits/stdc++.h> 
using namespace std;
#include <iostream>
#include <windows.h>
#include <cstdlib>
#include <fstream>
#include <filesystem>
#include <random>
#include <unordered_set>
#include <boost/random/mersenne_twister.hpp>
#include <boost/graph/random.hpp>

int main()
{
	string pathloop = "C:\\Users\\murat\\Desktop\\sensornetwork\\networkoptim\\datalist.txt";
	string line;
	
	boost::random::uniform_int_distribution<int> fdist(50000, 60000);
	boost::random::uniform_int_distribution<int> Ddist(100, 200);
	boost::random::uniform_int_distribution<int> cdist(1, 100);
	
	string upperfolder = "forest fire";
	int count = 40;
	ofstream seedtofile("C:\\Users\\murat\\Desktop\\sensornetwork\\networkoptim\\seedlist" + upperfolder+".txt");

	ifstream fileloop(pathloop);
	if (fileloop.is_open())
	{
		while (getline(fileloop, line))
		{
			int N;
			string filename;
			istringstream Streamloop(line);
			Streamloop >> filename >> N;

			boost::mt19937 gen(count);	
			seedtofile << filename << "\t" << count << "\n";

			//append cij
			string path = "C:\\Users\\murat\\Desktop\\sensornetwork\\networkoptim\\input\\"+upperfolder+"\\" + filename + ".txt";

			// Open "trial.txt" for reading
			ifstream inputFile(path);

			if (inputFile.is_open())
			{
				// Read the content of the file line by line
				std::vector<std::string> lines;
				std::string line;
				while (std::getline(inputFile, line))
					lines.push_back(line);

				inputFile.close();			

				std::ofstream openedfile(path);

				if (openedfile.is_open())
				{
					// Append text to each line and write the modified content back to the file
					for (const std::string& originalLine : lines)
					{
						std::string modifiedLine = originalLine + " " + to_string(cdist(gen));
						openedfile << modifiedLine << std::endl;
					}

					openedfile.close();
				}
			}

			std::string filePath = "C:\\Users\\murat\\Desktop\\sensornetwork\\networkoptim\\input\\" + upperfolder + "\\" + filename;

			// Convert string to dynamic char array
			std::vector<char> folderPath(filePath.begin(), filePath.end());
			folderPath.push_back('\0'); // Null-terminate the string

			if (CreateDirectoryA(folderPath.data(), NULL)) 
				std::cout << "Folder created successfully." << std::endl;
			else 
				std::cerr << "Failed to create the folder." << std::endl;

			ofstream ftofile("C:\\Users\\murat\\Desktop\\sensornetwork\\networkoptim\\input\\" + upperfolder + "\\" + filename + "\\f.txt");
			ofstream dtofile("C:\\Users\\murat\\Desktop\\sensornetwork\\networkoptim\\input\\" + upperfolder + "\\" + filename + "\\d.txt");

			int* f = new int[N];
			int* d = new int[N];

			for (int i = 0; i < N; ++i)
			{
				d[i] = Ddist(gen);
				dtofile << d[i] << "\n";
				f[i] = fdist(gen);
				ftofile << f[i] << "\n";
			}

			//MOVE THE txt file to newly created folder

			std::string sourceFile = path; // Replace with your source file name
			std::string destinationFolder = filePath; // Replace with your destination folder name

			if (CreateDirectoryA(destinationFolder.c_str(), NULL) || ERROR_ALREADY_EXISTS == GetLastError()) {
				std::string destinationPath = destinationFolder + "\\" + filename + ".txt";
				if (MoveFileA(sourceFile.c_str(), destinationPath.c_str())) {
					std::cout << "File moved successfully." << std::endl;
				}
				else {
					std::cerr << "Failed to move the file. Error code: " << GetLastError() << std::endl;
				}
			}
			else {
				std::cerr << "Failed to create the destination folder." << std::endl;
			}

			delete f;
			delete d;
			++count;
		}
	}
		
	return 0;
}