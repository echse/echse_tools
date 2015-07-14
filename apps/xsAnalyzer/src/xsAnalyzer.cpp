#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "system/system.h"
#include "except/except.h"
#include "cmdline/cmdline.h"
#include "table/table.h"
#include "typeconv/typeconv.h"
#include "mathutils/mathutils.h"
#include "xsections/xsections.h"

using namespace std;

int main (const int argc, const char* argv[]) {

  // Constants
  const string chars_colsep="\t ";
  const string chars_comment="#";
  const string colname_offset="offset";
  const string colname_x="x";
  const string colname_y="y";
  const string colname_elev="z";
  const string colname_rough="Kst";

  const double heightprecision= 0.001;
  const unsigned int maxiterations= 1000;

  // User input variables
  string file_xs;           // Input file with XS data
  string file_q;            // Input file with flow data
  string file_out;          // Output file
  bool read3d;              // Read input as 3d (x,y,z) data?
  bool readrough;           // Roughness data contained in file_xs?
  double defaultrough;      // Default roughness to use (inverse of Manning's n)
  double slope;             // Bed slope (-)
  double plain_length;      // Length of channel (m), projection in horizontal plane
  bool print_only_routing;  // Reduce output to data relevant for routing procedures?

  // Other variables
  xsection xs;
  vector<string> strings;
  vector<double> q, stage, topwidth, wetperim, wetarea, volume, dhdq, dvdq;
  vector<unsigned int> overtopped;
  double q_upper, q_lower, dhdq_upper, dhdq_lower, dvdq_upper, dvdq_lower;
  double characteristic_length;
  unsigned int nreserv;
  unsigned int nflows;
  ofstream ostr;

  try {

    ////////////////////////////////////////////////////////////////////////////
    // Get control data
    ////////////////////////////////////////////////////////////////////////////
    cout << "Getting command line input ..." << endl;
    cmdline cmdargs(argc, argv);
    try {
      file_xs= cmdargs.get_value("file_xsection", "=");
      file_q= cmdargs.get_value("file_flows", "=");
      file_out= cmdargs.get_value("file_out", "=");
      read3d= (cmdargs.get_value(
        "read_3d", "=").find_first_of("tTyY1") != string::npos);
      readrough= (cmdargs.get_value(
        "read_roughness", "=").find_first_of("tTyY1") != string::npos);
      defaultrough= as_double(cmdargs.get_value("default_roughness", "="));
      slope= as_double(cmdargs.get_value("slope", "="));
      plain_length= as_double(cmdargs.get_value("plain_length", "="));
      print_only_routing= (cmdargs.get_value(
        "print_only_routing", "=").find_first_of("tTyY1") != string::npos);
    } catch (except) {
      except e(__PRETTY_FUNCTION__,"Failed to retrieve info from command line.",
        __FILE__,__LINE__);
      throw(e);
    }
    // Check input data
    cout << "Checking input ..." << endl;
    if ((slope <= 0.) || (slope > 0.5)) {
      except e(__PRETTY_FUNCTION__,"Slope (given as dy/dx) not reasonable.",__FILE__,__LINE__);
      throw(e);
    }
    if (plain_length <= 0.) {
      except e(__PRETTY_FUNCTION__,"Reach length not reasonable.",__FILE__,__LINE__);
      throw(e);
    }
    if (!readrough && ((defaultrough <= 0.) || (defaultrough > 100.))) {
      except e(__PRETTY_FUNCTION__,"Used default roughness not reasonable.",__FILE__,__LINE__);
      throw(e);
    }

    ////////////////////////////////////////////////////////////////////////////
    // Read data
    ////////////////////////////////////////////////////////////////////////////
    cout << "Reading cross-section ..." << endl;
    try {
      xs.read(file_xs, chars_comment, chars_colsep,
        colname_offset, colname_x, colname_y, colname_elev, colname_rough, 
        "unknown", "unknown", -9999., read3d, readrough, defaultrough, -9999.);
    } catch (except) {
        stringstream errmsg;
        errmsg << "Could not read cross-section from file '" << file_xs << "'.";
        except e(__PRETTY_FUNCTION__,errmsg,__FILE__,__LINE__);
      throw(e);
    }
    cout << "Reading flow data ..." << endl;
    try {
      table qtab;
      qtab.read(file_q, false, chars_comment, chars_colsep);
      qtab.get_col(1, strings);
      nflows= qtab.nrow();
      q.resize(nflows);
      convert_type(strings, q); 
      qtab.clear();
    } catch (except) {
      stringstream errmsg;
      errmsg << "Could not read flow data from file '" << file_q << "'.";
      except e(__PRETTY_FUNCTION__,errmsg,__FILE__,__LINE__);
      throw(e);
    }
    if (nflows < 2) {
      stringstream errmsg;
      errmsg << "Number of flow values must be > 2. Found " << nflows << 
        " values in '" << file_q << "'.";
      except e(__PRETTY_FUNCTION__,errmsg,__FILE__,__LINE__);
      throw(e);
    }
    for (unsigned int i=0; i<(nflows-1); i++) {
      if (q[i+1] <= q[i]) {
        stringstream errmsg;
        errmsg << "Flow data in file '" << file_q << "' should be in ascending order.";
        except e(__PRETTY_FUNCTION__,errmsg,__FILE__,__LINE__);
        throw(e);
      }
    }

    ////////////////////////////////////////////////////////////////////////////
    // Compute simple xsection characteristics as a function of the flow
    ////////////////////////////////////////////////////////////////////////////
    cout << "Computing x-section properties ..." << endl;
    stage.resize(nflows);
    topwidth.resize(nflows);
    wetperim.resize(nflows);
    wetarea.resize(nflows);
    volume.resize(nflows);
    overtopped.resize(nflows);
    // Loop through flows
    for (unsigned int i=0; i<q.size(); i++) {
      // Get stage
      try {
        stage[i]= xs.return_surface(q[i], slope, heightprecision, maxiterations);
      } catch (except) {
        stringstream errmsg;
        errmsg << "Could not compute stage at normal depth for flow value " <<
          (i+1) << " of " << nflows << " (value=" << q[i] << ").";
        except e(__PRETTY_FUNCTION__,errmsg,__FILE__,__LINE__);
        throw(e);
      }
      // Check if overtopping occurs
      if ((stage[i] > xs.return_left_elevation()) |
          (stage[i] > xs.return_right_elevation())) {
        overtopped[i]= 1;
      } else {
        overtopped[i]= 0;        
      }

      // Get corresponding stage functions
      xs.return_stagefunctions(stage[i], wetarea[i], topwidth[i], wetperim[i]);
      // Compute volume
      volume[i]= wetarea[i] * plain_length;
      // Make sure that the values increase with rising flow
      if (i > 0) {
        if (stage[i] <= stage[i-1]) {
          stringstream errmsg;
          errmsg << "Stage does not increase for rising flow. For flow=" <<
            q[i-1] << " a stage of " << stage[i-1] << " was obtained but" <<
            " for flow=" << q[i] << " a stage of " << stage[i] << 
            " was computed. This may result from an improper/missing" <<
            " sub-division of the cross-section using a variable roughness.";
          except e(__PRETTY_FUNCTION__,errmsg,__FILE__,__LINE__);
          throw(e);
        }
      }
    }

    ////////////////////////////////////////////////////////////////////////////
    // Compute derivatives d_stage/d_flow and d_volume/d_flow
    ////////////////////////////////////////////////////////////////////////////
    dhdq.resize(nflows);
    dvdq.resize(nflows);

    // For flows > min and < max
    for (unsigned int i=1; i<(nflows-1); i++) {
      q_upper=    0.5*(q[i+1]+q[i]);
      dhdq_upper= (stage[i+1]-stage[i]) / (q[i+1]-q[i]);
      dvdq_upper= (volume[i+1]-volume[i]) / (q[i+1]-q[i]);
      q_lower=    0.5*(q[i]+q[i-1]);
      dhdq_lower= (stage[i]-stage[i-1]) / (q[i]-q[i-1]);
      dvdq_lower= (volume[i]-volume[i-1]) / (q[i]-q[i-1]);

      dhdq[i]= interpol(q_lower, q_upper, dhdq_lower, dhdq_upper, q[i]);
      dvdq[i]= interpol(q_lower, q_upper, dvdq_lower, dvdq_upper, q[i]);
    }
    // Values at margins
    dhdq[0]= (stage[1]-stage[0]) / (q[1]-q[0]);
    dvdq[0]= (volume[1]-volume[0]) / (q[1]-q[0]);
    dhdq[nflows-1]= (stage[nflows-1]-stage[nflows-2]) / (q[nflows-1]-q[nflows-2]);
    dvdq[nflows-1]= (volume[nflows-1]-volume[nflows-2]) / (q[nflows-1]-q[nflows-2]);

    ////////////////////////////////////////////////////////////////////////////
    // Write output
    ////////////////////////////////////////////////////////////////////////////
    cout << "Creating output ..." << endl;
    if (file_exists(file_out)) {
      stringstream errmsg;
      errmsg << "Output file '" << file_out << "' already exists.";
      except e(__PRETTY_FUNCTION__,errmsg,__FILE__,__LINE__);
      throw(e);
    }
    ostr.open(file_out.c_str());
    if (!ostr) {
      stringstream errmsg;
      errmsg << "Could not open output file '" << file_out << "'.";
      except e(__PRETTY_FUNCTION__,errmsg,__FILE__,__LINE__);
      throw(e);
    }
    ostr << "# Cross-section property table" << endl;
    ostr << "#   X-section data file:    " << file_xs << endl;
    ostr << "#   Flow data file:         " << file_q << endl;
    ostr << "#   Used reach length:      " << plain_length << endl;
    ostr << "#   Used slope:             " << slope << endl;
    if (!readrough) {
      ostr << "#   Used roughness: " << defaultrough << " (Global default)" << endl;
    } else {
      ostr << "#   Used roughness: Values supplied in x-section file" << endl;
    }
    if (!print_only_routing) {
      ostr << "flow" << "\t" << "stage" << "\t" << "overtop" << "\t" << "wet_area" << "\t" <<
        "top_width" << "\t" << "wet_perimeter" << "\t"  << "volume_total" <<
        "\t" << "dvdq_total" << "\t" << "nsub" << "\t"  << "volume_sub" <<
        "\t" << "dvdq_sub" << endl;
    } else {
      ostr << "flow" << "\t" << "overtop" << "\t" << "volume_total" << "\t" << "dvdq_total" <<
        "\t" << "nsub" << "\t"  << "volume_sub" << "\t" << "dvdq_sub" << endl;
    }
    for (unsigned int i=0; i<q.size(); i++) {
      // Compute characteristic length
      characteristic_length= q[i]/slope*dhdq[i];
      // Compute number of reservoirs
      nreserv= max(1, static_cast<int>(round(plain_length/characteristic_length)));
      // Write data
      if (!print_only_routing) {
        ostr <<
          numprint(q[i],3) << "\t" <<
          numprint(stage[i],3) << "\t" <<
          overtopped[i] << "\t" <<
          numprint(wetarea[i],2) << "\t" <<
          numprint(topwidth[i],2) << "\t" <<
          numprint(wetperim[i],2) << "\t" <<
          numprint(volume[i],3) << "\t" <<
          numprint(dvdq[i],3) << "\t" <<
          nreserv << "\t" <<
          numprint(volume[i]/nreserv,3) << "\t" <<
          numprint(dvdq[i]/nreserv,3) <<
          endl;
      } else {
        ostr <<
          numprint(q[i],3) << "\t" <<
          overtopped[i] << "\t" <<
          numprint(volume[i],3) << "\t" <<
          numprint(dvdq[i],3) << "\t" <<
          nreserv << "\t" <<
          numprint(volume[i]/nreserv,3) << "\t" <<
          numprint(dvdq[i]/nreserv,3) <<
          endl;
      }
    }
    ostr.close();

    ////////////////////////////////////////////////////////////////////////////
    // Clean up
    ////////////////////////////////////////////////////////////////////////////
    cout << "Cleaning up..." << endl;
    xs.clear();

    q.clear();
    stage.clear();
    topwidth.clear();
    wetperim.clear();
    wetarea.clear();
    volume.clear();
    dhdq.clear();
    dvdq.clear();
    overtopped.clear();


  // Terminate
	}
  catch (except) {
    except e(__PRETTY_FUNCTION__,"Stopped due to exception. See traceback.",__FILE__,__LINE__);
    e.print();
    return(1);
  }
  return(0);
}

