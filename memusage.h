#ifndef _MEMUSAGE_H_
#define _MEMUSAGE_H_

#include <iostream>
#include <sstream>
#include <string>
#include <fstream>


unsigned long get_memusage_vmpeak();
unsigned long get_memusage_vmsize();
void print_memusage(std::string pm_file_out);

#endif
