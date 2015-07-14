#include<iostream>
#include<fstream>
#include<sstream>
#include<string>

#include "typeconv/typeconv.h"
#include "table/table.h"
#include "cmdline/cmdline.h"
#include "except/except.h"
#include "arrays/arrays.h"

#include "meteofill_sub.h"

using namespace std;

int main(const int argc, const char* argv[]) {

  string ifile_data;
  string ifile_locations;
  double nodata;
  string ofile_data;
  string ofile_locations;
  string logfile;
  bool overwrite;
  string chars_colsep, chars_comment, char_colsep_out;
  double idw_power;
  unsigned int nsectors;
  unsigned int norigins;
  unsigned int ndigits_max;

  string line;
  unsigned int numline;
  string::size_type strpos;
  vector<string> tokens;
  vector<string> vect_id;
  vector<double> vect_x, vect_y, vect_z;
  table tab_locations;
  unsigned int col_id, col_x, col_y, col_z;
  bool ok;

  vector<string> vect_locations;
  vector<double> vect_values;
  vector<double> vect_values_old;
  vector<bool> is_missing;
  vector<bool> was_missing;
  bool config_changed;
  unsigned int numdata;

  vector<double> x_tar,y_tar,x_src,y_src;
  vector<unsigned int> inds_target,inds_source;
  vector<double> weights;
  double roundfac;

  // Variables for residual interpolation
  double limo_r2, limo_slope, limo_intercept;
  bool resid_apply;
  unsigned int resid_nmin;  // Min. number of data for trend estimation
  double resid_r2min; // Min. r^2 to apply residual interpolation
  double resid_llim;  // Lower limit of result of res. int. (avoids undershooting) 
  double resid_ulim;  // Upper limit of result of res. int. (avoids overshooting)


  try {

    // Get arguments
    cout << "Reading command line..." << endl;
    cmdline cmdargs(argc, argv);
    try {
      ifile_locations= cmdargs.get_value("ifile_locations", "=");
      ifile_data= cmdargs.get_value("ifile_data", "=");
      chars_colsep= cmdargs.get_value("chars_colsep", "=");
      chars_comment= cmdargs.get_value("chars_comment", "=");
      nodata= as_double(cmdargs.get_value("nodata", "="));
      idw_power= as_double(cmdargs.get_value("idw_power", "="));
      nsectors= as_unsigned_integer(cmdargs.get_value("nsectors", "="));
      norigins= as_unsigned_integer(cmdargs.get_value("norigins", "="));
      resid_nmin= as_unsigned_integer(cmdargs.get_value("resid_nmin", "="));
      resid_r2min= as_double(cmdargs.get_value("resid_r2min", "="));
      resid_llim= as_double(cmdargs.get_value("resid_llim", "="));
      resid_ulim= as_double(cmdargs.get_value("resid_ulim", "="));
      ofile_data= cmdargs.get_value("ofile_data", "=");
      ofile_locations= cmdargs.get_value("ofile_locations", "=");
      logfile= cmdargs.get_value("logfile", "=");
      ndigits_max= as_unsigned_integer(cmdargs.get_value("ndigits_max", "="));
      overwrite= as_logical(cmdargs.get_value("overwrite", "="));
    } catch (except) {
      except e(__PRETTY_FUNCTION__,"Failed to retrieve info from command line.",
        __FILE__,__LINE__);
      throw(e);
    }

    // Set column seperator for output
    char_colsep_out= chars_colsep.substr(0,1);

    // Set constant for rounding of interpolated data
    roundfac= pow(10.,ndigits_max);

    // Read table of locations
    cout << "Reading table of locations..." << endl;
    try {
      tab_locations.read(ifile_locations,true,chars_colsep,chars_comment);
    } catch (except) {
      stringstream errmsg;
      errmsg << "Failed to read table of locations from '" << ifile_locations << "'.";
      except e(__PRETTY_FUNCTION__,errmsg,__FILE__,__LINE__);
      throw(e);
    }
    try {
      col_id= tab_locations.colindex("id");
      col_x= tab_locations.colindex("x");
      col_y= tab_locations.colindex("y");
      col_z= tab_locations.colindex("z");
    } catch (except) {
      stringstream errmsg;
      errmsg << "Missing column(s) in file '" << ifile_locations << "'.";
      except e(__PRETTY_FUNCTION__,errmsg,__FILE__,__LINE__);
      throw(e);
    }

    // Open data output
    if (file_exists(ofile_data) && (!overwrite)) {
      stringstream errmsg;
      errmsg << "File '" << ofile_data << "' already exists.";
		  except e(__PRETTY_FUNCTION__,errmsg,__FILE__,__LINE__);
		  throw(e);
    }
    ofstream ofs_data;
    ofs_data.open(ofile_data.c_str());
    if (!ofs_data.is_open()) {
      stringstream errmsg;
      errmsg << "Unable to open file '" << ofile_data << "' for writing.";
		  except e(__PRETTY_FUNCTION__,errmsg,__FILE__,__LINE__);
		  throw(e);
    }

    // Open locations output
    if (file_exists(ofile_locations) && (!overwrite)) {
      stringstream errmsg;
      errmsg << "File '" << ofile_locations << "' already exists.";
		  except e(__PRETTY_FUNCTION__,errmsg,__FILE__,__LINE__);
		  throw(e);
    }
    ofstream ofs_locs;
    ofs_locs.open(ofile_locations.c_str());
    if (!ofs_locs.is_open()) {
      stringstream errmsg;
      errmsg << "Unable to open file '" << ofile_locations << "' for writing.";
		  except e(__PRETTY_FUNCTION__,errmsg,__FILE__,__LINE__);
		  throw(e);
    }
    ofs_locs << "id" << char_colsep_out << "x" << char_colsep_out << "y" <<
      char_colsep_out << "z" << endl;

    // Open logfile
    if (file_exists(logfile) && (!overwrite)) {
      stringstream errmsg;
      errmsg << "File '" << logfile << "' already exists.";
		  except e(__PRETTY_FUNCTION__,errmsg,__FILE__,__LINE__);
		  throw(e);
    }
    ofstream lg;
    lg.open(logfile.c_str());
    if (!lg.is_open()) {
      stringstream errmsg;
      errmsg << "Unable to open log file '" << logfile << "' for writing.";
		  except e(__PRETTY_FUNCTION__,errmsg,__FILE__,__LINE__);
		  throw(e);
    }
    lg << "time" << char_colsep_out << "numdata" << char_colsep_out <<
      "availability" << char_colsep_out << "config_changed" << char_colsep_out <<
      "residual_interp" << char_colsep_out << "r2" << endl;

    // Open data file
    cout << "Opening data file..." << endl;
    ifstream ifs;
    ifs.open(ifile_data.c_str());
    if (!ifs.is_open()) {
      stringstream errmsg;
      errmsg << "Cannot open file '" << ifile_data << "' for reading.";
      except e(__PRETTY_FUNCTION__,errmsg,__FILE__,__LINE__);
      throw(e);
    }
    numline= 0;

    // Try to read location IDs from first line of data file
    cout << "Reading location IDs from data file header..." << endl;
    while(!ifs.eof()) {
      getline(ifs,line);
      numline++;
	    strpos= line.find_first_not_of(" \t\n");    // Find first non-blank position
	    if (strpos != string::npos) {               // Skip blank/whitespace lines
        if (line.substr(strpos,1).find_first_of(chars_comment) == string::npos) { // Skip comment lines
          split(line, tokens, chars_colsep);
          break;
        }
	    }
    }
    if (tokens.size() < 2) {
      stringstream errmsg;
      errmsg << "Expecting at least two columns in file '" << ifile_data << "'.";
      except e(__PRETTY_FUNCTION__,errmsg,__FILE__,__LINE__);
      throw(e);
    }
    vect_locations.resize(tokens.size()-1);  // First column with times is skipped
    for (unsigned int i=1; i<tokens.size(); i++) {
      vect_locations[i-1]= tokens[i];
    }

    // Determine coordinates of the locations
    cout << "Assigning coordinates..." << endl;
    vect_x.resize(vect_locations.size());
    vect_y.resize(vect_locations.size());
    vect_z.resize(vect_locations.size());
    for (unsigned int i=0; i<vect_locations.size(); i++) {
      ok=false;
      for (unsigned int k=1; k<=tab_locations.nrow(); k++) {
        if (tab_locations.get_element(k,col_id) == vect_locations[i]) {
          try {
            vect_x[i]= as_double(tab_locations.get_element(k,col_x));
            vect_y[i]= as_double(tab_locations.get_element(k,col_y));
            vect_z[i]= as_double(tab_locations.get_element(k,col_z));
          } catch (except) {
            stringstream errmsg;
            errmsg << "Cannot extract coordinates (as double) from table" <<
              " read from file '" << ifile_locations << "'.";
            except e(__PRETTY_FUNCTION__,errmsg,__FILE__,__LINE__);
            throw(e);
          }
          ok=true;
          break;
        }
      }
      if (!ok) {
        stringstream errmsg;
        errmsg << "Missing record for location with ID '" << vect_locations[i] <<
          "' in coordinates table read from file '" << ifile_locations << "'.";
        except e(__PRETTY_FUNCTION__,errmsg,__FILE__,__LINE__);
        throw(e);
      }
      // Add location to output table
      ofs_locs << vect_locations[i] << char_colsep_out << numprint(vect_x[i],3) <<
        char_colsep_out << numprint(vect_y[i],3) << char_colsep_out <<
        numprint(vect_z[i],3) << endl;
    }

    // Print header
    ofs_data << tokens[0];
    for (unsigned int i=0; i<vect_locations.size(); i++) {
      ofs_data << char_colsep_out << vect_locations[i];
    }
    ofs_data << endl;

    // Read and process remaining lines
    cout << "Processing time series data..." << endl;
    // Init
    vect_values.resize(vect_locations.size());
    is_missing.resize(vect_locations.size());
    was_missing.resize(vect_locations.size());
    for (unsigned int i=0; i<was_missing.size(); i++) was_missing[i]= false;
    vect_values_old.resize(vect_locations.size());
    for (unsigned int i=0; i<vect_values_old.size(); i++) {
      vect_values_old[i]= nodata;
    }
    // Loop through records
    while(!ifs.eof()) {
      getline(ifs,line);
      numline++;
	    strpos= line.find_first_not_of(" \t\n");    // Find first non-blank position
	    if (strpos != string::npos) {               // Skip blank/whitespace lines
        if (line.substr(strpos,1).find_first_of(chars_comment) == string::npos) { // Skip comment lines
          // Check number of tokens
          split(line, tokens, chars_colsep);
          if (tokens.size() != (1+vect_locations.size())) {
            stringstream errmsg;
            errmsg << "Number of columns at line " << numline << " of file '" <<
              ifile_data << "' does not match the number of columns in" <<
              " the file header.";
            except e(__PRETTY_FUNCTION__,errmsg,__FILE__,__LINE__);
            throw(e);
          }
          // Convert values to double (first column with times is skipped)
          for (unsigned int i=1; i<tokens.size(); i++) {
            try {
              vect_values[i-1]= as_double(tokens[i]);
            } catch (except) {
              stringstream errmsg;
              errmsg << "Token at line " << numline << ", column " << i << " of file '" <<
                ifile_data << "' could not be converted into a numeric value.";
              except e(__PRETTY_FUNCTION__,errmsg,__FILE__,__LINE__);
              throw(e);
            }
          }
          // Has the data availability changed?
          numdata= 0;
          for (unsigned int i=0; i<vect_values.size(); i++) {
            if (vect_values[i] == nodata) {
              is_missing[i]= true;
            } else {
              is_missing[i]= false;
              numdata++;
            }
          }
          config_changed= false;
          for (unsigned int i=0; i<vect_values.size(); i++) {
            if (is_missing[i] != was_missing[i]) {
              config_changed= true;
              break;
            }
          }

          // Estimate missing data for this time step

          // Fill with old values, if there are no data at any location
          if (numdata == 0) {
            for (unsigned int i=0; i<vect_values.size(); i++) {
              if (vect_values_old[i] == nodata) {
                stringstream errmsg;
                errmsg << "Cannot substitute missing data at line " <<
                  numline << " of file '" << ifile_data << "' by older data." <<
                  " Older data are not available.";
                except e(__PRETTY_FUNCTION__,errmsg,__FILE__,__LINE__);
                throw(e);
              } else {
                vect_values[i]= vect_values_old[i];
              }
            }
            resid_apply= false;
  
          // Fill by spatial interpolation if data at some (but not all) locations exist
          } else if (numdata < vect_values.size()) {
            // Compute new weights if necessary
            if (config_changed) {
              station_weights(vect_x,vect_y,is_missing,idw_power,nsectors,norigins,
                inds_target,inds_source,weights);
            }

            // Compute linear trend using the stations' z-coordinate as a predictor
            if (numdata >= resid_nmin) {
              try {
                ok= linearModel(vect_z, vect_values, nodata,
                  limo_slope, limo_intercept, limo_r2);
                resid_apply= ok & (limo_r2 >= resid_r2min);
              } catch (except) {
                stringstream errmsg;
                errmsg << "Failed to compute linear regression.";
                except e(__PRETTY_FUNCTION__,errmsg,__FILE__,__LINE__);
                throw(e);
              }
            } else {
              resid_apply= false;
            }

            // Apply residual interpolation
            if (resid_apply) {
              // Initialize missing data using the linear trend
              for (unsigned int i=0; i<is_missing.size(); i++) {
                if (is_missing[i]) vect_values[i]= vect_z[i] * limo_slope + limo_intercept;
              }
              // Add weighted residuals
              for (unsigned int i=0; i<inds_target.size(); i++) {
                vect_values[inds_target[i]]= vect_values[inds_target[i]] + 
                  (vect_values[inds_source[i]] -
                  (vect_z[inds_source[i]] * limo_slope + limo_intercept)) * weights[i];
              }
              // Correct undesired extrapolation
              for (unsigned int i=0; i<inds_target.size(); i++) {
                if (vect_values[inds_target[i]] < resid_llim)
                  vect_values[inds_target[i]]= resid_llim;
                if (vect_values[inds_target[i]] > resid_ulim)
                  vect_values[inds_target[i]]= resid_ulim;
              }
            // Normal interpolation of original values
            } else {
              // Initialize missing data to zero
              for (unsigned int i=0; i<is_missing.size(); i++) {
                if (is_missing[i]) vect_values[i]= 0.;
              }
              // Add weighted values
              for (unsigned int i=0; i<inds_target.size(); i++) {
                vect_values[inds_target[i]]= vect_values[inds_target[i]] + 
                  vect_values[inds_source[i]] * weights[i];
              }
            }

            // Round interpolated values
            for (unsigned int i=0; i<is_missing.size(); i++) {
              if (is_missing[i]) vect_values[i]= round(vect_values[i] * roundfac) / roundfac;
            }
            
          } else {
            // else: Data complete, nothing to do
            resid_apply= false;
          }
    
          // Save values of this step for possible extrapolation in time
          for (unsigned int i=0; i<vect_values.size(); i++) {
            vect_values_old[i]= vect_values[i];
          }

          // Save missing locations to detect changes in station config
          for (unsigned int i=0; i<is_missing.size(); i++) {
            was_missing[i]= is_missing[i];
          }
        
          // Add record to output
          ofs_data << tokens[0];
          for (unsigned int i=0; i<vect_values.size(); i++) {
            ofs_data << char_colsep_out << vect_values[i];
          }
          ofs_data << endl;

          // Add log info
          lg << tokens[0] << char_colsep_out << numdata << char_colsep_out <<
            (numdata * 100. / vect_values.size()) << char_colsep_out <<
            config_changed << char_colsep_out;
          if (resid_apply) {
            lg << "true" << char_colsep_out << limo_r2 << endl;
          } else {
            lg << "false" << char_colsep_out << "NA" << endl;
          }

        }
	    }
    } // End of loop over time steps

    // Close output files
    ofs_data.close();
    ofs_locs.close();
    lg.close();

    // Clean up
    tab_locations.clear();
    vect_locations.clear();
    vect_x.clear();
    vect_y.clear();
    vect_z.clear();
    
    cout << "Done." << endl;


  } catch (except) {
    except e(__PRETTY_FUNCTION__,"Stopped due to error.",__FILE__,__LINE__);
    e.print();
    return(1);
  }
  return(0);
}

