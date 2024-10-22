#include <fstream>   // for file operations
#include <iostream>  // for cout

void loadFromFile(std::string filename)
{
  // open file
  std::ifstream file(filename.c_str(), std::ios::in);

  // check if file is open
  if (!file.is_open())
  {
    std::cout << "Could not open parameter file \"" << filename << "\"." << std::endl;
    return;
  }

  // loop over lines of file
  for (int lineNo = 0;; lineNo++)
  {
    // read line
    std::string line;
    getline(file, line);

    // at the end of the file break for loop
    if (file.eof())
      break;

    // print line
    std::cout << "line " << lineNo << ": " << line << std::endl;
  }
}