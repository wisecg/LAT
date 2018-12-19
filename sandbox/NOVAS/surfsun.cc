// Compile:
// g++ -o GenerateSolarPosition `root-config --cflags` novas.c novascon.c nutation.c readeph0.c solsys3.c surfsun.cc `root-config --glibs`
// This code uses NOVAS and various correction factors to calculate the angle of the sun
// Usage: Just run the damn thing!
// If you don't understand time: http://stjarnhimlen.se/comp/time.html
// Accuracy of everything is calculated on a per-minute basis, code can be modified to be more accurate

#include <stdio.h>
#include <stdlib.h>
#include "novas.h"
#include <iostream>
#include <fstream>
#include "TFile.h"
#include "TTree.h"

using namespace std;

int main (int argc, char** argv)
{
   int Year = 0, Month = 0, Day = 0;
   double Hour=0, Minute=0;
   double UTC[60], Azimuth[60], Zenith[60], SolarPos[60];

   // Difference between UT1 and UTC, based off of rotation of the Earth
   // Obtained from http://maia.usno.navy.mil/search/search.html, changes every day
   // If the search website is down, use the ftp ftp://maia.usno.navy.mil/ser7/
   double ut1_utc = 0;


   TFile *f1 = new TFile("SunTreeUTC_Full.root", "RECREATE");
   TTree *tree = new TTree("SunTree", "Solar Angles calculated wrt SURF");
   tree->Branch("Year", &Year, "Year/I");
   tree->Branch("Month", &Month, "Month/I");
   tree->Branch("Day", &Day, "Day/I");
   tree->Branch("Hour", &Hour, "Hour/D");
   tree->Branch("ut1_utc", &ut1_utc, "ut1_utc/D");
   tree->Branch("UTC", &UTC, "UTC[60]/D");
   tree->Branch("Azimuth", &Azimuth, "Azimuth[60]/D");
   tree->Branch("Zenith", &Zenith, "Zenith[60]/D");
   tree->Branch("SolarPos", &SolarPos, "SolarPos[60]/D");

   short int error = 0;
   short int accuracy = 1;
   // coord_sys (short int)
   //    Code specifying coordinate system of the output position.
   //       = 0 ... GCRS or "local GCRS"
   //       = 1 ... true equator and equinox of date
   //       = 2 ... true equator and CIO of date
   //       = 3 ... astrometric coordinates, i.e., without light
   //               deflection or aberration.
	short int coord_sys = 1;

   // Observer
   on_surface geo_loc;
   observer obs_loc;

   const double latitude = 44.353;
   // const double longitude = 103.767;
   const double longitude = -103.767; // Correct longitude is actually negative
   const double height = 110.642; // Has no effect
   const double temperature = 10.0; // Has no effect
   const double pressure = 1010.0; // Has no effect

   make_on_surface(latitude, longitude, height, temperature, pressure, &geo_loc);
   make_observer_on_surface(latitude, longitude, height, temperature, pressure, &obs_loc);

   // Make sun -- first create dummy catalog entry, then create object
   cat_entry dummy_star;
   object sun;
   make_cat_entry("DUMMY"," ",0,0.0,0.0,0.0,0.0,0.0,0.0, &dummy_star);
   if ((error = make_object(0, 10, "Sun", &dummy_star, &sun)) != 0)
   {
      printf ("Error %d from make_object (Sun)\n", error);
      return (error);
   }

   const double timezone = -0.0;

   // Correction terms
   // Leap seconds obtained from http://maia.usno.navy.mil/ser7/tai-utc.dat, changes every year
   // If the search website is down, use the ftp ftp://maia.usno.navy.mil/ser7/
   short int leap_secs = 0;

   sky_pos t_place;
   double jd_utc, jd_tt, jd_ut1, jd_tdb, delta_t = 0;
   double rat, dect, dist;
   double zd = 0, az = 0, rar = 0, decr = 0;

   ifstream input;
   input.open("UT1_UTC.txt");

   int dYear, dMonth, dDay;
   double dCorr;

   int dLoop = 0;
   while(!input.eof())
   {
      dLoop++;
      if(dLoop%100==0) cout << "Loop: " << dLoop << endl;
      input >> dYear;
      input >> dMonth;
      input >> dDay;
      input >> dCorr;

      Year = 2000 + dYear;
      Month = dMonth;
      Day = dDay;
      ut1_utc = dCorr;

      if(Year >= 2012 && Year < 2015) leap_secs = 35;
      else if(Year == 2015 || Year == 2016) leap_secs = 36;
      else if(Year == 2017 || Year == 2018) leap_secs = 37;
      // 2019 and 2020 should have 38 leap seconds
      else {
         cout << "Year not in database, using 38 seconds as correction" << endl;
         leap_secs = 38;
      }

      // 32.184 corrects for terrestrial time (tt)
      delta_t = 32.184 + leap_secs - ut1_utc;

      for(int i = 0; i < 24; i++)
      {
         Hour = i;

         for(int j = 0; j < 60; j++)
         {
            jd_utc = julian_date(Year, Month, Day, i + ((double)j/60.) + timezone);
            jd_tt = jd_utc + ((double) leap_secs + 32.184) / 86400.0;
            jd_ut1 = jd_utc + ut1_utc / 86400.0;
            jd_tdb = jd_tt;          /* Approximation good to 0.0017 seconds. */

            if ((error = place(jd_tt, &sun, &obs_loc, delta_t, coord_sys, accuracy, &t_place)) != 0)
            {
               printf ("Error %d from place.", error);
               return (error);
            }

            equ2hor (jd_ut1, delta_t, accuracy, 0.0, 0.0, &geo_loc, t_place.ra , t_place.dec, 1, &zd, &az, &rar, &decr);
            UTC[j] = jd_utc;
            Zenith[j] = zd;
            Azimuth[j] = az;
            SolarPos[j] = t_place.dis;
         }

         tree->Fill();
      }


   }

   f1->cd();
   tree->Write();
   f1->Close();

   return 0;
}
