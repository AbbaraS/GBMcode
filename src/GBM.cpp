// GBM.cpp : Defines the entry point for the application.
//

#include "GBM.h"

namespace GBM
{
	std::string GBM::ReadColumnHeaders(const std::string& filename)
	{
		std::ifstream file(filename);
		if (!file.is_open())
		{
			std::cerr << "Could not open file: " << filename << std::endl;
			exit(0);
			return "";
		}
		if (file.good())
		{
			std::string line;
			std::getline(file, line);
			return line;
		}
		else
		{
			std::cerr << "Could not read file: " << filename << std::endl;
			exit(0);
			return "";
		}
	}

	std::vector<std::string> GBM::ReadUserInput(const std::string& filename)
	{
		std::ifstream file(filename);
		if (!file.is_open())
		{
			std::cerr << "Could not open file: " << filename << std::endl;
			exit(0);
			return {};
		}
		std::vector<std::string> features;
		std::string line;
		while (std::getline(file, line))
		{
			if (!line.empty())
			{
				// Remove newlines
				line.erase(std::remove(line.begin(), line.end(), '\r'), line.cend());
				line.erase(std::remove(line.begin(), line.end(), '\n'), line.cend());
				// Remove any quotes from the line
				line.erase(std::remove(line.begin(), line.end(), '\"'), line.cend());
			}
			features.push_back(line);
		}
		return features;
	}

	std::vector<std::string> GBM::ExtractPatientIDs(const std::string& filename)
	{
		std::ifstream file(filename);
		if (!file.is_open())
		{
			std::cerr << "Could not open file: " << filename << std::endl;
			exit(0);
			return {};
		}
		auto geneIds = GBM::ReadColumnHeaders(filename);
		std::vector<std::string> patientIds;
		size_t pos = 0;

		// The first three sections of the gene ID are the patient ID
		std::regex pattern(R"((\w+-\d+-\d+))");
		std::smatch match;

		std::string geneId;
		// The column headers are tab separated
		while ((pos = geneIds.find("\t")) != std::string::npos)
		{
			// Extract the gene ID
			geneId = geneIds.substr(0, pos);
			if (std::regex_search(geneId, match, pattern)) {
				patientIds.push_back(match.str());
			}
			geneIds.erase(0, pos + 1);
		}
		// Search last bit of the string
		if (std::regex_search(geneIds, match, pattern)) {
			patientIds.push_back(match.str());
		}
		return patientIds;
	}

	Table GBM::ReadFileAsTable(const std::string& filename)
	{
		std::ifstream file(filename);
		if (!file.is_open())
		{
			std::cerr << "Could not open file: " << filename << std::endl;
			exit(0);
			return {};
		}
		Table table{};

		std::string line;
		bool firstLine = true;
		// Get row
		while (std::getline(file, line))
		{
			line.erase(std::remove(line.begin(), line.end(), '\r'), line.cend());
			line.erase(std::remove(line.begin(), line.end(), '\n'), line.cend());
			line.erase(std::remove(line.begin(), line.end(), '\"'), line.cend());

			std::stringstream linestream(line);
			std::string cell;

			Row row{};

			if (firstLine)
			{
				row.cells.push_back(""); // Empty cell for the header
				firstLine = false;
			}
			// Get cell
			while (std::getline(linestream, cell, ','))
			{
				row.cells.push_back(cell);
			}
			table.AddRow(row);
		}
		return table;
	}

	int GBM::Table::GetHeaderIndex(const std::string& header)
	{
		for (size_t i = 0; i < rows[0].cells.size(); i++)
		{
			if (rows[0].cells[i] == header)
			{
				return i;
			}
		}
		return -1;
	}

	std::vector<int> GBM::Table::GetPatientRow(const std::string& patientId)
	{
		// Get index of patient barcode
		auto barcodeIndex = GetHeaderIndex("bcr_patient_barcode");
		std::vector<int> indecies;
		for (size_t i = 0; i < rows.size(); i++)
		{
			if (rows[i].cells[barcodeIndex] == patientId)
			{
				indecies.push_back(i);
			}
		}
		return indecies;
	}

	std::vector<std::string> GBM::Table::GetCells(const std::string& header, const std::string& patientId)
	{
		auto headerIndex = GetHeaderIndex(header);
		auto patientIndex = GetPatientRow(patientId);
		if (headerIndex == -1 || patientIndex.size() == 0)
		{
			return {};
		}
		std::vector<std::string> cells;
		for (auto& index : patientIndex)
		{
			cells.push_back(rows[index].cells[headerIndex]);
		}
		return cells;
	}

	bool GBM::Table::WriteToFile(const std::string& filename)
	{
		std::ofstream file(filename);
		if (!file.is_open())
		{
			std::cerr << "Could not open file: " << filename << std::endl;
			exit(0);
			return false;
		}
		for (auto& row : rows)
		{
			for (size_t i = 0; i < row.cells.size(); i++)
			{
				file << row.cells[i];
				if (i != row.cells.size() - 1)
				{
					file << ",";
				}
			}
			file << std::endl;
		}
		return true;
	}
}

int main()
{
	using namespace GBM;
	// Read in patient IDs and all tables
	auto patientIds = ExtractPatientIDs("GBM_RNAseqdata_HTSEQ_FKPM.harmonized.txt");
	auto drugTable = ReadFileAsTable("clinical_drug_GBM.txt");
	auto patientTable = ReadFileAsTable("clinical_patient_GBM.txt");
	auto followupTable = ReadFileAsTable("clinical_followup_GBM.txt");

	std::vector<Table> tables{ drugTable, patientTable, followupTable };

	// Read in user input
	std::cout << "Enter the location of the user input file: ";
	std::string userInputLocation;
	std::cin >> userInputLocation;
	auto userInput = ReadUserInput(userInputLocation);
	std::cout << "Enter the location to save the output file: ";
	std::string outputLocation;
	std::cin >> outputLocation;
	
	Table outTable{};
	Row headers{};
	headers.cells.push_back("Patient ID");
	// Push the features to the headers
	headers.cells.insert(headers.cells.end(), userInput.begin(), userInput.end());
	outTable.AddRow(headers);

	for (auto& patient : patientIds)
	{
		Row row;
		// Just get the four digit number for Patient ID
		row.cells.push_back(patient.substr(patient.find_last_of('-') + 1));
		bool hasValue = false;
		for (auto& header : headers.cells)
		{
			if (header == "Patient ID")
			{
				continue;
			}
			std::vector<std::string> allCells;
			for (auto& table : tables)
			{
				auto cells = table.GetCells(header, patient);
				if (cells.size() == 0)
				{
					continue;
				}
				allCells.insert(allCells.end(), cells.begin(), cells.end());
			}
			if (allCells.size() == 0)
			{
				row.cells.push_back("");
			}
			else
			{
				std::string fullCell;
				for (size_t i = 0; i < allCells.size(); i++)
				{
					fullCell += allCells[i];
					if (i != allCells.size() - 1)
					{
						fullCell += " ";
					}
				}
				row.cells.push_back(fullCell);
				hasValue = true;
			}
		}
		if (hasValue)
		{
			outTable.AddRow(row);
		}
	}
	outTable.WriteToFile(outputLocation);
	std::cout << "Output file saved to: " << outputLocation << std::endl;
	return 0;
}
