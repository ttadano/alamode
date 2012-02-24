#include "parsephon.h"
#include <iostream>
#include <string>

using namespace PHON_NS;

Input::Input(PHON *phon, int narg, char **arg): Pointers(phon) {}

Input::~Input() {}

void Input::parce_input()
{
    using namespace std;
    string job_title, mode;

    int dos, non_analytic;

    cin >> job_title;
}