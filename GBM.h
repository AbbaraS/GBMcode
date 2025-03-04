// suleima_test.h : Include file for standard system include files,
// or project specific include files.

#pragma once

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <regex>
#include <sstream>

namespace GBM
{
	struct Row
	{
		std::vector<std::string> cells;
	};

	class Table
	{
	public:
		Table() = default;

		void			 AddRow(Row row) { rows.push_back(row); }
		int				 GetHeaderIndex(const std::string& header);
		std::vector<int> GetPatientRow(const std::string& patientId);
		std::vector<std::string>		 GetCells(const std::string& header, const std::string& patientId);
		bool			 WriteToFile(const std::string& filename);

	private:
		std::vector<Row> rows;
	};

	Table ReadFileAsTable(const std::string& filename);
	std::string ReadColumnHeaders(const std::string& filename);

	std::vector<std::string> ReadUserInput(const std::string& filename);

	std::vector<std::string> ExtractPatientIDs(const std::string& filename);
}

// TODO: Reference additional headers your program requires here.
