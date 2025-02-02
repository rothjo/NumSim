#include "output_writer_serial/output_writer_text.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

void OutputWriterText::writeFile(double currentTime)
{
  // Assemble the filename
  std::stringstream fileName;
  fileName << "out/output_" << std::setw(4) << std::setfill('0') << fileNo_ << ".txt";
  
  // increment file no.
  fileNo_++;

  // open file
  std::ofstream file(fileName.str().c_str());

  if (!file.is_open())
  {
    std::cout << "Could not write to file \"" << fileName.str() << "\".";
    return;
  }

  // write time
  file << "t: " << currentTime << std::endl;

  // write mesh width
  file << "nCells: " << discretization_->nCells()[0] << "x" << discretization_->nCells()[1] 
    << ", dx: " << discretization_->dx() << ", dy: " << discretization_->dy() << std::endl << std::endl;

  const int fieldWidth = 9;   // number of characters to use for a single value

  // write u
  // ---------
  // write header lines
  file << "u (" << discretization_->u().size()[0] << "x" << discretization_->u().size()[1] << "): " << std::endl 
    << std::string(fieldWidth, ' ') << "|";
  for (int i = discretization_->uIBegin() - 1; i < discretization_->uIEnd() + 1; i++)
  {
    file << std::setw(fieldWidth) << i;
  }
  file << std::endl << std::string(fieldWidth*(discretization_->u().size()[0]+2)+1, '-') << std::endl;

  // write u values
  for (int j = discretization_->uJEnd(); j >= discretization_->uJBegin() - 1; j--)
  {
    file << std::setw(fieldWidth) << j << "|";
    for (int i = discretization_->uIBegin() - 1; i < discretization_->uIEnd() + 1; i++)
    {
      file << std::setw(fieldWidth) << std::setprecision(fieldWidth-6) << discretization_->u(i,j);
    }
    file << std::endl;
  }
  file << std::endl;

  // write v
  // ---------
  // write header lines
  file << "v (" << discretization_->v().size()[0] << "x" << discretization_->v().size()[1] << "): " << std::endl 
    << std::string(fieldWidth, ' ') << "|";
  for (int i = discretization_->vIBegin() - 1; i < discretization_->vIEnd() + 1; i++)
  {
    file << std::setw(fieldWidth) << i;
  }
  file << std::endl << std::string(fieldWidth*(discretization_->v().size()[0]+2)+1, '-') << std::endl;

  // write v values
  for (int j = discretization_->vJEnd() + 1; j >= discretization_->vJBegin() - 1; j--)
  {
    file << std::setw(fieldWidth) << j << "|";
    for (int i = discretization_->vIBegin() - 1; i < discretization_->vIEnd() + 1; i++)
    {
      file << std::setw(fieldWidth) << std::setprecision(fieldWidth-6) << discretization_->v(i,j);
    }
    file << std::endl;
  }
  file << std::endl;

  // write p
  // ---------
  // write header lines
  file << "p (" << discretization_->p().size()[0] << "x" << discretization_->p().size()[1] << "): " << std::endl 
    << std::string(fieldWidth, ' ') << "|";
  for (int i = discretization_->pIBegin() - 1; i < discretization_->pIEnd()+1; i++)
  {
    file << std::setw(fieldWidth) << i;
  }
  file << std::endl << std::string(fieldWidth*(discretization_->p().size()[0]+2)+1, '-') << std::endl;

  // write p values
  for (int j = discretization_->pJEnd(); j >= discretization_->pJBegin() - 1; j--)
  {
    file << std::setw(fieldWidth) << j << "|";
    for (int i = discretization_->pIBegin() - 1; i < discretization_->pIEnd() + 1; i++)
    {
      file << std::setw(fieldWidth) << std::setprecision(fieldWidth-6) << discretization_->p(i,j);
    }
    file << std::endl;
  }
  file << std::endl;

  // write f
  // ---------
  // write header lines
  file << "F (" << discretization_->u().size()[0] << "x" << discretization_->u().size()[1] << "): " << std::endl 
    << std::string(fieldWidth, ' ') << "|";
  for (int i = discretization_->uIBegin() - 1; i < discretization_->uIEnd() + 1; i++)
  {
    file << std::setw(fieldWidth) << i;
  }
  file << std::endl << std::string(fieldWidth*(discretization_->u().size()[0]+2)+1, '-') << std::endl;

  // write f values
  for (int j = discretization_->uJEnd(); j >= discretization_->uJBegin() - 1; j--)
  {
    file << std::setw(fieldWidth) << j << "|";
    for (int i = discretization_->uIBegin() - 1; i < discretization_->uIEnd() + 1; i++)
    {
      file << std::setw(fieldWidth) << std::setprecision(fieldWidth-6) << discretization_->f(i,j);
    }
    file << std::endl;
  }
  file << std::endl;

  // write g
  // ---------
  // write header lines
  file << "G (" << discretization_->v().size()[0] << "x" << discretization_->v().size()[1] << "): " << std::endl 
    << std::string(fieldWidth, ' ') << "|";
  for (int i = discretization_->vIBegin() - 1; i < discretization_->vIEnd() + 1; i++)
  {
    file << std::setw(fieldWidth) << i;
  }
  file << std::endl << std::string(fieldWidth*(discretization_->v().size()[0]+2)+1, '-') << std::endl;

  // write g values
  for (int j = discretization_->vJEnd(); j >= discretization_->vJBegin() - 1; j--)
  {
    file << std::setw(fieldWidth) << j << "|";
    for (int i = discretization_->vIBegin() - 1; i < discretization_->vIEnd() + 1; i++)
    {
      file << std::setw(fieldWidth) << std::setprecision(fieldWidth-6) << discretization_->g(i,j);
    }
    file << std::endl;
  }
  file << std::endl;

  // write rhs
  // ---------
  // write header lines
  file << "rhs (" << discretization_->p().size()[0] << "x" << discretization_->p().size()[1] << "): " << std::endl 
    << std::string(fieldWidth, ' ') << "|";
  for (int i = discretization_->pIBegin() - 1; i < discretization_->pIEnd() + 1; i++)
  {
    file << std::setw(fieldWidth) << i;
  }
  file << std::endl << std::string(fieldWidth*(discretization_->p().size()[0]+2)+1, '-') << std::endl;

  // write rhs values
  for (int j = discretization_->pJEnd(); j >= discretization_->pJBegin()-1; j--)
  {
    file << std::setw(fieldWidth) << j << "|";
    for (int i = discretization_->pIBegin()-1; i < discretization_->pIEnd()+1; i++)
    {
      file << std::setw(fieldWidth) << std::setprecision(fieldWidth-6) << discretization_->rhs(i,j);
    }
    file << std::endl;
  }
  file << std::endl;

}

void OutputWriterText::writePressureFile()
{
  // counter for files, counter value is part of the file name
  static int pressurefileNo = 0;

  // Assemble the filename
  std::stringstream fileName;
  fileName << "out/pressure_" << std::setw(4) << std::setfill('0') << pressurefileNo++ << ".txt";
  
  // open file
  std::ofstream file(fileName.str().c_str());

  if (!file.is_open())
  {
    std::cout << "Could not write to file \"" << fileName.str() << "\".";
    return;
  }

  // write mesh width
  file << "nCells: " << discretization_->nCells()[0] << "x" << discretization_->nCells()[1] 
    << ", dx: " << discretization_->dx() << ", dy: " << discretization_->dy() << std::endl << std::endl;

  const int fieldWidth = 9;   // number of characters to use for a single value

  // write p
  // ---------
  // write header lines
  file << "p (" << discretization_->p().size()[0] << "x" << discretization_->p().size()[1] << "): " << std::endl 
    << std::string(fieldWidth, ' ') << "|";
  for (int i = discretization_->pIBegin(); i < discretization_->pIEnd(); i++)
  {
    file << std::setw(fieldWidth) << i;
  }
  file << std::endl << std::string(fieldWidth*(discretization_->p().size()[0]+2)+1, '-') << std::endl;

  // write p values
  for (int j = discretization_->pJEnd()-1; j >= discretization_->pJBegin(); j--)
  {
    file << std::setw(fieldWidth) << j << "|";
    for (int i = discretization_->pIBegin(); i < discretization_->pIEnd(); i++)
    {
      file << std::setw(fieldWidth) << std::setprecision(fieldWidth-6) << discretization_->p(i,j);
    }
    file << std::endl;
  }
  file << std::endl;

}
