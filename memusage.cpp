#include <cstring>
#include <unistd.h>
#include <iostream>
#include "memusage.h"


using namespace std;

unsigned long get_memusage_vmpeak()
{
    pid_t curr_pid = getpid();
    stringstream  ss_fn_in;
    ss_fn_in << "/proc/" << curr_pid << "/status";
    string  str_fn_in;
    ss_fn_in >> str_fn_in;
    ifstream in(str_fn_in.c_str());

    if (!in.is_open()){
            cerr << "ERROR: cannot open file " << str_fn_in << endl;
            exit(1);
    }

    char  line [4096];

    char  vname[128];
    unsigned long vsize;

    do {
        in.getline(line, 4095);
        sscanf(line, "%s %lu", vname, &vsize);
    } while (strcmp(vname, "VmPeak:") != 0);

    return vsize;
}

unsigned long get_memusage_vmsize()
{
    pid_t curr_pid = getpid();
    stringstream  ss_fn_in;
    ss_fn_in << "/proc/" << curr_pid << "/status";
    string  str_fn_in;
    ss_fn_in >> str_fn_in;
    ifstream in(str_fn_in.c_str());

    if (!in.is_open()){
            cerr << "ERROR: cannot open file " << str_fn_in << endl;
            exit(1);
    }

    char  line [4096];

    char  vname[128];
    unsigned long vsize;

    do {
        in.getline(line, 4095);
        sscanf(line, "%s %lu", vname, &vsize);
    } while (strcmp(vname, "VmSize:") != 0);

    return vsize;
}

void print_memusage(string pm_file_out)
{
    pid_t curr_pid = getpid();
    stringstream  ss_fn_in;
    ss_fn_in << "/proc/" << curr_pid << "/status";
    string  str_fn_in;
    ss_fn_in >> str_fn_in;
    ifstream in(str_fn_in.c_str());

    if (!in.is_open()){
            cerr << "ERROR: cannot open file " << str_fn_in << endl;
            exit(1);
    }

    char  line [4096];

    char  vname[128];
    unsigned long vsize;

    do {
        in.getline(line, 4095);
        sscanf(line, "%s %lu", vname, &vsize);
        cout << "vname = " << vname << " vsize = " << vsize << endl;
    } while (strcmp(vname, "VmPeak:") != 0);
    cout << line << endl;

    int tw;
    cout << "Type an int to continue:" << endl;
    cin >> tw;
}
