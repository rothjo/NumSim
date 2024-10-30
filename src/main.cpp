#include "output_writer/write_paraview_output.h"
#include "settings.h"

#include <iostream>
#include <cstdlib>

int main(int argc, char *argv[])
{
  // write 5 output files
  // for (int i = 0; i < 5; i++)
  // {
  //   writeParaviewOutput(i);
  // }
  Settings settings;
  settings.loadFromFile("parameters.txt");
  settings.printSettings();

  return EXIT_SUCCESS;
}