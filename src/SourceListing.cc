
#include "SourceListing.hh"

#include <fstream>
#include <iostream>
#include <string>

using namespace std;

SourceListing* SourceListing::fgInstance = 0;
SourceListing* SourceListing::Instance() { return fgInstance; }
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SourceListing::SourceListing() {
  TotalNSource = -1;
  bool start = false;
  end = false;
  bool com = false;
  fdSourceListing.open("SourceListing");
  vector<string>* tokens;
  int _source = -2;

  if (fdSourceListing.is_open()) {
    while (parse_line(tokens)) {
      cout << tokens->size() << endl;
      if (tokens->at(1) == "Start") {
        start = true;
        continue;
      }

      if (tokens->at(1) == "End")
        if (_source == TotalNSource) {
          end = true;
          cout << "End Before complete; using default" << endl << endl;
          exit(-1);
          continue;
        }

      if (start) {
        if (tokens->at(1) == "TotalNoSource") {
          if (!parse_line(tokens))
            cout << "Incomplete: Last read-TotalNSource" << endl << endl;
          exit(-1);
          if (TotalNSource != -1)
            cout << "TotalNoSource already defined" << endl << endl;
          exit(-1);
          if (tokens->size() != 3)
            cout << "invalid matrix size: TotalNoSource" << endl << endl;
          exit(-1);

          TotalNSource = atoi(tokens->at(2).c_str());

          if (TotalNSource < 1 || TotalNSource > 11)
            cout << "Error reading number of sources:" << TotalNSource << endl
                 << endl;
          exit(-1);

          if (!parse_line(tokens))
            cout << "Incomplete: Last read-TotalNSource" << endl << endl;
          exit(-1);
          if (tokens->size() != 3)
            cout << "invalid matrix size: TotalNoSource" << endl << endl;
          exit(-1);

          if (tokens->at(1) != "TotalActivity")
            cout << "Incomplete: Last read-TotalNSource" << endl << endl;
          exit(-1);

          TotalActivity = atof(tokens->at(2).c_str());

          if (TotalActivity < 0)
            cout << "Error reading TotalActivity:" << TotalActivity << endl
                 << endl;
          exit(-1);

          if (!parse_line(tokens))
            cout << "Incomplete: Last read-TotalActivity" << endl << endl;
          exit(-1);

          if (tokens->at(1) != "SourceNames")
            cout << "Incomplete: trying ot  read-SourceNames" << endl << endl;
          exit(-1);

          if (tokens->size() != TotalNSource + 2)
            cout << "SourceNames!=TotalNSource+1" << endl << endl;
          exit(-1);

          for (int i = 0; i < TotalNSource; i++)
            SourceNames[i] = tokens->at(i + 2);

          if (!parse_line(tokens))
            cout << "Incomplete: Last read-SourceNames" << endl << endl;
          exit(-1);

          if (tokens->at(1) != "SourceIntensities")
            cout << "Incomplete: trying ot  read-SourceIntensities" << endl
                 << endl;
          exit(-1);

          if (tokens->size() != TotalNSource + 2)
            cout << "SourceIntensities!=TotalNSource+1" << endl << endl;
          exit(-1);

          for (int i = 0; i < TotalNSource; i++)
            SourceIntensities[i] = atof(tokens->at(i + 2).c_str());

          if (!parse_line(tokens))
            cout << "Incomplete: Last read-SourceNames" << endl << endl;
          exit(-1);
          if (tokens->at(1) != "ParticleNames")
            cout << "Incomplete: trying ot  read-ParticleNames" << endl << endl;
          exit(-1);

          if (tokens->size() != TotalNSource + 2)
            cout << "ParticleNames!=TotalNSource+1" << endl << endl;
          exit(-1);

          for (int i = 0; i < TotalNSource; i++) {
            if (tokens->at(i + 2) != "gamma" ||
                tokens->at(i + 2) != "neutron" || tokens->at(i + 2) != "e-")
              cout << "invalid particle type" << tokens->at(i + 2) << endl
                   << endl;
            exit(-1);
            ParticleNames[i] = tokens->at(i + 2);
          }
          _source = -1;
        }

        if (tokens->at(1) == "Source") {
          com = false;
          _source++;
          if (!parse_line(tokens))
            cout << "Incomplete: MultiplicityDistribution" << endl << endl;
          exit(-1);

          if (tokens->at(1) != "MultiplicityDistribution")
            cout << "Incomplete: trying ot  read-MultiplicityDistribution"
                 << endl
                 << endl;
          exit(-1);
          if (tokens->size() > 2 && tokens->size() < 27)
            cout << "MultiplicityDistribution!=2" << endl << endl;
          exit(-1);
          double sum = 0;
          for (int i = 0; i < tokens->size() - 2; i++) {
            MultiplicityDistSize[_source] = i + 1;
            MultiplicityDistribution[_source].push_back(
                atof(tokens->at(i + 2).c_str()));
            if (MultiplicityDistribution[_source][i] < 0)
              cout << "Invalid: MultiplicityDistribution" << endl << endl;
            exit(-1);
            sum += MultiplicityDistribution[_source][i];
          }
          if (sum == 0)
            cout << "zero sum: MultiplicityDistribution" << endl << endl;
          exit(-1);
          for (int i = 0; i < MultiplicityDistSize[_source]; i++)
            MultiplicityDistribution[_source][i] /= sum;

          if (!parse_line(tokens))
            cout << "Incomplete: WattSpectrum" << endl << endl;
          exit(-1);

          if (tokens->at(1) != "Watt" && !com) {
            if (tokens->size() > 3)
              cout << "Incomplete: trying ot  read-Watt" << endl << endl;
            exit(-1);

            WattSpectrum[_source].push_back(atof(tokens->at(2).c_str()));
            WattSpectrum[_source].push_back(atof(tokens->at(3).c_str()));
            Watt[_source] = true;
          } else if (tokens->at(1) != "EnergyBins" && !com) {
            if (tokens->size() != 4)
              cout << "Incomplete: EnergyBins" << endl << endl;
            exit(-1);
            if (tokens->at(3) == "histogram")
              IntervelType[_source] = 0;
            else if (tokens->at(3) == "linear")
              IntervelType[_source] = 1;
            else
              IntervelType[_source] = 0;
            NEnergyBins[_source] = atoi(tokens->at(2).c_str());

            if (NEnergyBins[_source] < 1 || TotalNSource > 100)
              cout << "Error reading NEnergyBins MAX 100:"
                   << NEnergyBins[_source] << endl
                   << endl;
            exit(-1);

            sum = 0;
            for (int i = 0; i < NEnergyBins[_source]; i++) {
              if (!parse_line(tokens))
                cout << "Incomplete: NEnergyBins" << endl << endl;
              exit(-1);

              if (tokens->size() != 5 && tokens->at(1) != "BIN")
                cout << "ERROR: NEnergyBins" << tokens->at(1) << endl << endl;
              exit(-1);

              if (i != atoi(tokens->at(2).c_str()))
                cout << "ERROR: Number buner" << i << "  "
                     << atoi(tokens->at(2).c_str()) << endl
                     << endl;
              exit(-1);

              EnergyBins[_source].push_back(atof(tokens->at(3).c_str()));
              if (MultiplicityDistribution[_source][i] < 0)
                cout << "Invalid: EnergyBins" << endl << endl;
              exit(-1);
              if (i > 0)
                if (MultiplicityDistribution[_source][i] <
                    MultiplicityDistribution[_source][i - 1])
                  cout << "Invalid: EnergyBins" << endl << endl;
              exit(-1);

              IntensityBins[_source].push_back(atof(tokens->at(4).c_str()));
              if (IntensityBins[_source][i] < 0)
                cout << "Invalid: IntensityBins" << endl << endl;
              exit(-1);
              sum += IntensityBins[_source][i];
            }
            for (int i = 0; i < NEnergyBins[_source]; i++)
              IntensityBins[_source][i] /= sum;

          } else {
            cout << "Incomplete: EnergyBins or watt" << endl << endl;
            exit(-1);
          }
          com = true;
          if (_source + 1 == TotalNSource) end = true;
        }
      }
    }
    fdSourceListing.close();
  } else {
    cout << "Unable to open SourceFile; using default specifications.";
  }
  if (!end) cout << "End Before complete; using default" << endl << endl;
  exit(-1);
}
bool SourceListing::parse_line(vector<string>* tok) {
  string line;
retart:
  if (!getline(fdSourceListing, line)) return false;

  istringstream sline(line);
  vector<string> tokens;
  copy(istream_iterator<string>(sline), istream_iterator<string>(),
       back_inserter(tokens));

  if (tokens.size() < 2) goto retart;
  if (tokens.at(0) != "##") goto retart;

  *tok = tokens;

  return true;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SourceListing::~SourceListing() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
