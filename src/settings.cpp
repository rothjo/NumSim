#include "settings.h"
#include <fstream>
#include <iomanip>

/**
 * Parse a text file with settings, each line contains "<parameterName> = <value>"
 * 
 * @param filename name of the file to parse
 * @return void
 */
void Settings::loadFromFile(std::string filename)
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
    // check if end of file is reached
    if (file.eof()) {
      break;
    }
    // read line
    std::string line;
    getline(file, line);

    // remove whitespace at beginning of line
    line.erase(0, line.find_first_not_of(" \t")); 

    // if first character is a '#', skip line
    if (line[0] == '#')
      continue; 

    // if line does not contain a '=' sign, skip line
    if (line.find('=') == std::string::npos)
      continue; 

    // extract and isolate parameter name
    std::string parameterName = line.substr(0, line.find('='));
    if (parameterName.find_first_of(" \t") != std::string::npos)
    {
      parameterName.erase(parameterName.find_first_of(" \t"));
    }

    // extract and isolate value
    std::string value = line.substr(line.find('=') + 1);

    // remove whitespace at beginning of value
    value.erase(0, value.find_first_not_of(" \t"));
    // remove comments at end of value
    if (value.find_first_of("#") != std::string::npos)
      value.erase(value.find_first_of("#"));
    // remove whitespace at end of value
    value.erase(value.find_last_not_of(" \t") + 1); 
    
    // match parameterName and value to given variables
    if (parameterName == "endTime"){
      endTime = atof(value.c_str());
    }
    else if (parameterName == "physicalSizeX"){ 
      physicalSize[0] = atof(value.c_str());
    }
    else if (parameterName == "physicalSizeY"){
      physicalSize[1] = atof(value.c_str());
    }
    else if (parameterName == "re"){
      re = atof(value.c_str());
    }
    else if (parameterName == "gX"){
      g[0] = atof(value.c_str());
    }
    else if (parameterName == "gY"){
      g[1] = atof(value.c_str());
    }
    else if (parameterName == "dirichletBottomX"){
      dirichletBcBottom[0] = atof(value.c_str());
    }
    else if (parameterName == "dirichletBottomY"){
      dirichletBcBottom[1] = atof(value.c_str());
    }
    else if (parameterName == "dirichletTopX"){
      dirichletBcTop[0] = atof(value.c_str());
    }
    else if (parameterName == "dirichletTopY"){
      dirichletBcTop[1] = atof(value.c_str());
    }
    else if (parameterName == "dirichletLeftX"){
      dirichletBcLeft[0] = atof(value.c_str());
    }
    else if (parameterName == "dirichletLeftY"){
      dirichletBcLeft[1] = atof(value.c_str());
    }
    else if (parameterName == "dirichletRightX"){
      dirichletBcRight[0] = atof(value.c_str());
    }
    else if (parameterName == "dirichletRightY"){
      dirichletBcRight[1] = atof(value.c_str());
    }
    else if (parameterName == "nCellsX"){
      nCells[0] = atoi(value.c_str());
    }
    else if (parameterName == "nCellsY"){
      nCells[1] = atoi(value.c_str());
    }
    else if (parameterName == "useDonorCell"){
      if(value == "true"){
        useDonorCell = true;
      }
      else if(value == "false"){
        useDonorCell = false;
      }
      else{
        std::cout << "Unknown value for useDonorCell: " << value << std::endl;
      }
    }
    else if (parameterName == "alpha"){
      alpha = atof(value.c_str());
    }
    else if (parameterName == "tau"){
      tau = atof(value.c_str());
    }
    else if (parameterName == "maximumDt"){
      maximumDt = atof(value.c_str());
    }
    else if (parameterName == "pressureSolver"){
      pressureSolver = value;
    }
    else if (parameterName == "omega"){
      omega = atof(value.c_str());
    }
    else if (parameterName == "epsilon"){
      epsilon = atof(value.c_str());
    }
    else if (parameterName == "maximumNumberOfIterations"){
      maximumNumberOfIterations = static_cast<int>(atof(value.c_str()));
    }
    else{
      std::cout << "Unknown parameter \"" << parameterName << "\" in line " << lineNo << std::endl;
    }
  }
}

/**
 * Output all settings to console
 * 
 * @return void
 */
void Settings::printSettings()
{
  std::cout << "Settings: " << std::endl
    << "  physicalSize: " << physicalSize[0] << " x " << physicalSize[1] << ", nCells: " << nCells[0] << " x " << nCells[1] << std::endl
    << "  endTime: " << endTime << " s, re: " << re << ", g: (" << g[0] << "," << g[1] << "), tau: " << tau << ", maximum dt: " << maximumDt << std::endl
    << "  dirichletBC: bottom: (" << dirichletBcBottom[0] << "," << dirichletBcBottom[1]  << ")"
    << ", top: ("  << dirichletBcTop[0] << "," << dirichletBcTop[1]  << ")"
    << ", left: ("  << dirichletBcLeft[0] << "," << dirichletBcLeft[1] << ")"
    << ", right: ("  << dirichletBcRight[0] << "," << dirichletBcRight[1] << ")" << std::endl
    << "  useDonorCell: " << std::boolalpha << useDonorCell << ", alpha: " << alpha << std::endl
    << "  pressureSolver: " << pressureSolver << ", omega: " << omega << ", epsilon: " << epsilon << ", maximumNumberOfIterations: " << maximumNumberOfIterations << std::endl;
}