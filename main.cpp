#include <iostream>
#include <cstdlib>

int main(int argc, char *argv[])
{
  // if the number of given command line arguments is only 1 (= the program name), print out usage information and exit
  if (argc == 1)
  {
    std::cout << "usage: " << argv[0] << " <filename>" << std::endl;

    return EXIT_FAILURE;
  }

  // read in the first argument
  std::string filename = argv[1];

  // print message
  std::cout << "Filename: \"" << filename << "\"" << std::endl;

  return EXIT_SUCCESS;
}