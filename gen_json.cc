#include <iostream>
#include <fstream>
#include <ctime>
#include "GATDataSet.hh"
#include "DataSetInfo.hh"
using namespace std;

void createJson(bool blind);

int main(int argc, const char** argv)
{
  bool blind=0;

  // parse user args
  vector<string> opt(argv + 1, argv + argc);
  for (size_t i = 0; i < opt.size(); i++) {
    if (opt[i] == "-b") { blind=1; cout << "Generating blind run sets...\n"; }
    if (opt[i] == "-h") {
      cout << "Usage: ./genJSON [options] : generate a JSON file for DataSetInfo run ranges.\n"
           << "       [-b : generate JSON file for blind data subsets]\n"
           << "       [-h : print this message]\n";
    }
  }

  // generate json file from DataSetInfo.cc
  createJson(blind);
}

void createJson(bool blind)
{
  time_t now = time(NULL);
  struct tm * curtime = localtime(&now);
  string fDate = asctime(curtime);
  fDate.erase(remove(fDate.begin(), fDate.end(), '\n'), fDate.end()); // remove the newline

  string fName = "runsBkg.json";
  if (blind) fName = "runsBlind.json";
  ofstream outFile(fName);

  vector<int> dataSets = {0,1,2,3,4,5,6};
  if (blind) {
    dataSets = {1,2,5,6};
  }
  cout << "Creating JSON run range file for Data Sets: ";
  for (auto i : dataSets) cout << i << " ";
  cout << endl;

  outFile << "{\n";
  outFile << Form("\"note1\":\"Generated from DataSetInfo.cc on %s\",\n",fDate.c_str());

  for (size_t i = 0; i < dataSets.size(); i++)
  {
    int dsNum = dataSets[i];
    GATDataSet ds;
    map<int, vector<int>> ranges;
    if (blind) LoadBlindDataSet(ds, dsNum, -1, ranges);
    else LoadDataSet(ds, dsNum, -1, ranges);

    // convert the map to JSON format
    outFile << Form("\"%i\": {\n", dsNum);
    size_t idx = 0;
    for (auto &j : ranges)
    {
      string sub = Form("  \"%i\": [",j.first);
      for (auto &k : j.second) {
        sub += Form("%i, ", k);
      }
      sub.erase(sub.length()-2);

      if (idx==ranges.size()-1)
        sub += "]\n";
      else
        sub += "],\n";

      outFile << sub;
      idx++;
    }
    if (i == dataSets.size()-1)
      outFile << "}\n";
    else
      outFile << "},\n";


  }
  outFile << "}\n";
  outFile.close();
}