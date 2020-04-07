//START REGION RELEASE
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "psrsalsa.h"
#include <cpgplot.h>

//START REGION DEVELOP

#define PgplotLibMaxStringLength 1000

//START REGION RELEASE

/* Variables set by pgplotMap and used by pgplotMapCoordinate */
static double internal_pgplot_xmin = 0;
static double internal_pgplot_xmax = 1;
static double internal_pgplot_ymin = 0;
static double internal_pgplot_ymax = 1;
static int internal_pgplot_nrx = 1;
static int internal_pgplot_nry = 1;

void print_pgplot_version_used(FILE *stream)
{
  char version[25];
  int length = 20;
  cpgqinf("VERSION", version, &length);
  fprintf(stream, "%s (library)", version);
}


/* Function sets nx and ny to the array ellements in the data that has
   been plotted with pgplotMap at coordinates (x, y). If a point is
   selected outside the array the point is set at the border. In that
   case the return value is 1 instead of 0. */
int pgplotMapCoordinate_dbl(double x, double y, int *nx, int *ny)
{
  int clip;
  clip = 0;
  *nx = 0.5 + (x - internal_pgplot_xmin)*(internal_pgplot_nrx-1.0)/(internal_pgplot_xmax-internal_pgplot_xmin);
  *ny = 0.5 + (y - internal_pgplot_ymin)*(internal_pgplot_nry-1.0)/(internal_pgplot_ymax-internal_pgplot_ymin);
  if(*nx < 0) {
    *nx = 0;
    clip = 1;
  }
  if(*nx >= internal_pgplot_nrx) {
    *nx = internal_pgplot_nrx - 1;
    clip = 1;
  }
  if(*ny < 0) {
    *ny = 0;
    clip = 1;
  }
  if(*ny >= internal_pgplot_nry) {
    *ny = internal_pgplot_nry - 1;
    clip = 1;
  }
  return clip;
}

int pgplotMapCoordinate(float x, float y, int *nx, int *ny)
{
  double x2, y2;
  x2 = x;
  y2 = y;
  return pgplotMapCoordinate_dbl(x2, y2, nx, ny);
}

//START REGION DEVELOP
//START REGION RELEASE

/* Inverse of pgplotMapCoordinate. Returns coordinates of the centre of bin. */
void pgplotMapCoordinateInverse_dbl(double *x, double *y, int nx, int ny)
{
  if(internal_pgplot_nrx != 1)
    *x = ((double)nx)*(internal_pgplot_xmax-internal_pgplot_xmin)/((double)internal_pgplot_nrx-1.0) + internal_pgplot_xmin;
  else
    *x = internal_pgplot_xmin;
  if(internal_pgplot_nry != 1)
    *y = ((double)ny)*(internal_pgplot_ymax-internal_pgplot_ymin)/((double)internal_pgplot_nry-1.0) + internal_pgplot_ymin;
  else
    *y = internal_pgplot_ymin;
}

void pgplotMapCoordinateInverse(float *x, float *y, int nx, int ny)
{
  double x2, y2;
  pgplotMapCoordinateInverse_dbl(&x2, &y2, nx, ny);
  *x = x2;
  *y = y2;
}

//START REGION DEVELOP
//START REGION RELEASE

/* Returns the size of one bin in the map plotted with pgplotMap. */
void pgplotMapCoordinateBinSize(float *dx, float *dy)
{
  *dx = (internal_pgplot_xmax-internal_pgplot_xmin)/((double)internal_pgplot_nrx-1.0);
  *dy = (internal_pgplot_ymax-internal_pgplot_ymin)/((double)internal_pgplot_nry-1.0);
}

//START REGION DEVELOP
//START REGION RELEASE

/* windowwidth & windowheight = (if specified) sets the resolution of the device. 
otherwise aspectratio (if set) sets the aspect ratio of the device, without specifying the resolution. */
void pgplot_setWindowsize(int windowwidth, int windowheight, float aspectratio)
{
  float x, y;
  if(windowwidth > 0 && windowheight > 0) {
    x = windowwidth*0.01175548589341692789994673739445429916373;
    y = (windowheight-1)/(float)windowwidth;
    ppgpap(x, y);
  }else if(aspectratio > 0) {
    ppgpap(0, aspectratio);
  }
}

//START REGION DEVELOP
//START REGION RELEASE

/* Sets the viewport and the axes ranges. It will also calculate TR,
   which relies on the internal_pgplot_XXXX parameters to be set. */
void pgplot_makeframe(pgplot_frame_def_internal *frame)
{
  if(frame->svp)
    ppgsvp(frame->svp_x1, frame->svp_x2, frame->svp_y1, frame->svp_y2);
  if(frame->swin) {
    ppgswin(frame->swin_x1, frame->swin_x2, frame->swin_y1, frame->swin_y2);
    //    printf("YYYYY %e %e\n", internal_pgplot_xmin,internal_pgplot_xmax);
    if(internal_pgplot_nrx != 1) {
      frame->TR[0] = internal_pgplot_xmin - (internal_pgplot_xmax-internal_pgplot_xmin)/(float)(internal_pgplot_nrx-1);  
      frame->TR[1] = (internal_pgplot_xmax-internal_pgplot_xmin)/(float)(internal_pgplot_nrx-1);
    }else {
      if(internal_pgplot_xmin != internal_pgplot_xmax) {
	frame->TR[0] = internal_pgplot_xmin - (internal_pgplot_xmax-internal_pgplot_xmin); 
	frame->TR[1] = internal_pgplot_xmax-internal_pgplot_xmin;
      }else {
	frame->TR[0] = internal_pgplot_xmin - (frame->swin_x2 - frame->swin_x1);  // Assume the bin width is equal to the range plotted. Note that since internal_pgplot_xmin = internal_pgplot_xmax, we don't have the information about what the binwidth is supposed to be.
	frame->TR[1] = frame->swin_x2 - frame->swin_x1;
      }
    }
    frame->TR[2] = 0;
    frame->TR[3] = internal_pgplot_ymin - (internal_pgplot_ymax-internal_pgplot_ymin)/(float)(internal_pgplot_nry-1); frame->TR[4] = 0;    frame->TR[5] = (internal_pgplot_ymax-internal_pgplot_ymin)/(float)(internal_pgplot_nry-1);
    if(internal_pgplot_nry == 1) {
      frame->TR[3] = internal_pgplot_ymin; 
      frame->TR[5] = (internal_pgplot_ymax-internal_pgplot_ymin);
      if(frame->TR[5] == 0) {             /* A hack for pplot if only 1 vector. There is something dodgy about all this code */
	frame->TR[3] = -1;
	frame->TR[5] = 1;
      }
    }
  }
}

//START REGION DEVELOP
//START REGION RELEASE

/*
Make a PA plot
data = data file
showtotpol = nonzero: plot the total degree of polarization
nopaswing = nonzero: do not produce PA-swing panel
showEll = nonzero: plot the ellipticity angle
viewport = The name of the device etc
xlabel, ylabel, ylabel_pa = labels (set title in pgplotbox (keywords are substituted with header parameters with the str_replace_header_params() function), xlabel and ylabel are ignored in pgplotbox)
longitude_left & longitude_right = pulse longitude range in degrees
xunit_type: 0=degrees, 1=phase
Imin, Imax = yrange of intensities. Set both to zero for autoscale
pa_bottom & pa_top = PA range
PAoffset = extra vertical offset PA's in plot
sigma_limit = required S/N to plot PA's
datalinewidth = line width of data
ysize2 = relative ysize of PA plot [1=default]
dashed = dashed lines instead (usefull for B/W plots)
noynumbers = Don't plot numbers on vertical axis
textoption = the option in command line for the texts, i.e. -text
herrorbaroption = the option in command line for the horizontal errorbars, i.e. -herrorbar
herrorbaroption2 = the option in command line for the horizontal errorbars in the pa plot, i.e. -herrorbar2
verrorbaroption = the option in command line for the vertical errorbars, i.e. -verrorbar
verrorbaroption2 = the option in command line for the vertical errorbars in the pa plot, i.e. -verrorbar2
argc, argv = number of command line options and the options themselves. 
             It parses for command lines like: 
             -text "x y ch lw font color" "text" 
outline_txt = The text appears as an outline
outline_lw, outline_color = line width and color of the outline
overlayPA = Overplot an RVM model
overlayalpha, overlaybeta, overlaypa0, overlayl0 = parameters of RVM curve
overlayPAfine: If set, a much finer resolution than the resolution of the data is used to draw the RVM curve
nrJumps: number of OPM jumps to be included in the RVM curve
jump_longitudes & jump_offsets: the list of longitudes and offsets where the OPMs occur
dontclose = plotdevice is not closed
dontopen = it assumes a plotting device is opened and selected
padist = the pa distribution to be attached underneath the plot, if not NULL pointer and if the data is already read in successfully
padist_pamin, padist_pamax = range in PA covered by PA distribution in degrees
padist_saturize: If > 1, the count rate will get saturated
elldist = same for ellipticity angle distribution
elldist_saturize: If > 1, the count rate will get saturated
nowedge: if non-zero, no wedge is shown next to PA/Ellipticity distributions
  Return 1: Succesful
         0: Plot device couldn't be opened or data couldn't be plotted
*/

int pgplotPAplot(datafile_definition data, int showtotpol, int nopaswing, int showEll, pgplot_options_definition *pgplot, char *xlabel, char *ylabel, char *ylabel_pa, char *ylabel_ell, float longitude_left, float longitude_right, int xunit_type, float Imin, float Imax, float pa_bottom, float pa_top, float PAoffset, float sigma_limit, float datalinewidth, float ysize2, int dashed, int noynumbers, char *textoption, char *herrorbaroption, char *herrorbaroption2, char *verrorbaroption, char *verrorbaroption2, int argc, char **argv, int outline_txt, int outline_lw, int outline_color, int overlayPA, float overlayalpha, float overlaybeta, float overlaypa0, float overlayl0, int overlayPAfine, int nrJumps, float *jump_longitudes, float *jump_offsets, datafile_definition *padist, float padist_pamin, float padist_pamax, float padist_saturize, datafile_definition *elldist, float elldist_saturize, int nowedge, verbose_definition verbose)
{
  int deviceID, text_ci, text_lw, text_f, ok, domove, showPAdist, showELLdist;
  float ymin, ymax, text_x, text_y, text_ch, I, Iold, overallplotscaling;
  long i, j;
  char *newtext;
  pgplot_frame_def_internal frame;
  pgplot_options_definition *pgplot_backup;
 
  int showwedge = 1;
  if(nowedge)
    showwedge = 0;

  showPAdist = 0;
//START REGION DEVELOP
  if(padist != NULL) {
    if(padist->format != MEMORY_format) {
      printerror(verbose.debug, "ERROR pgplotPAplot: PA distribution does not appear to be loaded into memory\n");
      return 0;      
    }
    if(padist->NrPols != 1) {
      printerror(verbose.debug, "ERROR pgplotPAplot: PA distribution is expected to have a single polarization channel\n");
      return 0;      
    }
    if(padist->gentype != GENTYPE_PADIST) {
      printwarning(verbose.debug, "WARNING pgplotPAplot: The file opened as a PA distribution does not have a gentype set to be a PA distribution. It is currently set to %s.", returnGenType_str(padist->gentype));
      return 0;      
    }
    if(padist->NrBins != data.NrBins) {
      printerror(verbose.debug, "ERROR pgplotPAplot: The number of pulse longitude bins is different in the PA distribution compared to the profile.");
      return 0;      
    }
    //    printf("XXXXX %f %f\n", padist_pamin, padist_pamax);
    if(padist_pamax - padist_pamin < 0) {
      float junk;
      junk = padist_pamax;
      padist_pamax = padist_pamin;
      padist_pamin = junk;
    }
    //    printf("XXXXX %f %f\n", padist_pamin, padist_pamax);
    if(padist_pamax - padist_pamin > 360) {
      printerror(verbose.debug, "ERROR pgplotPAplot: The PA distribution cannot cover more than 360 degrees.");
      return 0;      
    }
    showPAdist = 1;
  }
//START REGION RELEASE

  showELLdist = 0;
//START REGION DEVELOP
  if(elldist != NULL) {
    if(elldist->format != MEMORY_format) {
      printerror(verbose.debug, "ERROR pgplotPAplot: Ellipticity angle distribution does not appear to be loaded into memory\n");
      return 0;      
    }
    if(elldist->NrPols != 1) {
      printerror(verbose.debug, "ERROR pgplotPAplot: Ellipticity angle distribution is expected to have a single polarization channel\n");
      return 0;      
    }
    if(elldist->gentype != GENTYPE_ELLDIST) {
      printwarning(verbose.debug, "WARNING pgplotPAplot: The file opened as an ellipticity angle distribution does not have a gentype set to be an ellipticity angle distribution. It is currently set to %s.", returnGenType_str(elldist->gentype));
      return 0;      
    }
    if(elldist->NrBins != data.NrBins) {
      printerror(verbose.debug, "ERROR pgplotPAplot: The number of pulse longitude bins is different in the ellipticity angle distribution compared to the profile.");
      return 0;      
    }
    showELLdist = 1;
  }
//START REGION RELEASE

  overallplotscaling = 1;
  if(showEll && (showPAdist == 0 && showELLdist == 0)) { // Require another panel, so make all graphs a bit smaller
    overallplotscaling = 0.78;
  }else if(showEll && (showPAdist || showELLdist)) {
    overallplotscaling = 0.53;
  }else if(showEll == 0 && (showPAdist || showELLdist)) {
    overallplotscaling = 0.63;
  }
  if(showPAdist && showELLdist) {
    overallplotscaling = 0.4;
  }

  pgplot_backup = (pgplot_options_definition *)malloc(sizeof(pgplot_options_definition));
  if(pgplot_backup == NULL) {
    printerror(verbose.debug, "ERROR pgplotPAplot: Memory allocation error\n");
    return 0;      
  }
  memcpy(pgplot_backup, pgplot, sizeof(pgplot_options_definition));

  pgplot->box.box_labelsize *= overallplotscaling;
  pgplot->box.label_ch *= overallplotscaling;

  float vp_left, vp_right, vp_profile_top, vp_profile_bottom, vp_paswing_top, vp_paswing_bottom, vp_ellswing_top, vp_ellswing_bottom, vp_padist_top, vp_padist_bottom, vp_elldist_top, vp_elldist_bottom;
  vp_left            = 0.1+pgplot->viewport.dxplot;
  vp_right           = 0.1+pgplot->viewport.dxplot+0.8*pgplot->viewport.xsize*overallplotscaling;
  vp_profile_bottom  = 0.35+(1-overallplotscaling)*0.6+pgplot->viewport.dyplot;
  vp_profile_top     = vp_profile_bottom +0.55*pgplot->viewport.ysize*overallplotscaling;
  if(nopaswing == 0) {
    vp_paswing_bottom  = vp_profile_bottom      -0.25*pgplot->viewport.ysize*ysize2*overallplotscaling;
    vp_paswing_top     = vp_profile_top   - 0.55*pgplot->viewport.ysize*overallplotscaling;
  }else {
    vp_paswing_bottom  = vp_profile_bottom;
    vp_paswing_top     = vp_profile_bottom + overallplotscaling*0.25*pgplot->viewport.ysize*ysize2;
  }
  if(showEll) {
    vp_ellswing_bottom = vp_paswing_bottom  - overallplotscaling*0.25*pgplot->viewport.ysize*ysize2;
    vp_ellswing_top    = vp_paswing_top     - overallplotscaling*0.25*pgplot->viewport.ysize*ysize2;
  }else {
    vp_ellswing_bottom = vp_paswing_bottom;
    vp_ellswing_top    = vp_paswing_top;
  }
  if(showPAdist) {
    vp_padist_top    = vp_ellswing_bottom;
    vp_padist_bottom = vp_padist_top - overallplotscaling*0.5*pgplot->viewport.ysize*ysize2;
  }else {
    vp_padist_top    = vp_ellswing_top;
    vp_padist_bottom = vp_ellswing_bottom;
  }
  if(showELLdist) {
    vp_elldist_top    = vp_padist_bottom;
    vp_elldist_bottom = vp_elldist_top - overallplotscaling*0.5*pgplot->viewport.ysize*ysize2;
  }else {
    vp_elldist_top    = vp_padist_top;
    vp_elldist_bottom = vp_padist_bottom;
  }


  if(pgplot_opendevice(&(pgplot->viewport), &deviceID, verbose) == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR pgplotPAplot: Cannot open plot device");
    memcpy(pgplot, pgplot_backup, sizeof(pgplot_options_definition));
    free(pgplot_backup);
    return 0;
  }

  if(data.poltype == POLTYPE_PAdPA) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR pgplotPAplot: Data does only have PA points, but no profile");
    memcpy(pgplot, pgplot_backup, sizeof(pgplot_options_definition));
    free(pgplot_backup);
    return 0;
  }
  if(data.poltype != POLTYPE_ILVPAdPA && data.poltype != POLTYPE_ILVPAdPATEldEl) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR pgplotPAplot: Data does not appear to contain PA data");
    memcpy(pgplot, pgplot_backup, sizeof(pgplot_options_definition));
    free(pgplot_backup);
    return 0;
  }

  if(data.poltype == POLTYPE_ILVPAdPA && data.NrPols != 5) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR pgplotPAplot: Data does not appear to contain PA data as NrPols != 5");
    memcpy(pgplot, pgplot_backup, sizeof(pgplot_options_definition));
    free(pgplot_backup);
    return 0;
  }
  if(data.poltype == POLTYPE_ILVPAdPATEldEl && data.NrPols != 8) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR pgplotPAplot: Data does not appear to contain PA data as NrPols != 8");
    memcpy(pgplot, pgplot_backup, sizeof(pgplot_options_definition));
    free(pgplot_backup);
    return 0;
  }

  if((showtotpol || showEll) && data.poltype != POLTYPE_ILVPAdPATEldEl) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR pgplotPAplot: Data does not appear to contain the total polarization and/or ellipticity angle, but it was requested to be plotted");
    memcpy(pgplot, pgplot_backup, sizeof(pgplot_options_definition));
    free(pgplot_backup);
    return 0;
  }

  //Substitute header parameter in title keywords
  newtext = str_replace_header_params(data, pgplot->box.title, verbose);
  if(newtext == NULL) {
    printerror(verbose.debug, "ERROR pgplotPAplot: Cannot substuture header parameter in the tile");
    memcpy(pgplot, pgplot_backup, sizeof(pgplot_options_definition));
    free(pgplot_backup);
    return 0;
  }
  strcpy(pgplot->box.title, newtext);
  free(newtext);

  //  clear_pgplot_box(&box);
  clear_pgplot_frame(&frame);


  // Determine vertical range in top panel
  ymin = ymax = NAN;
  for(j = 1; j < (data.NrBins); j++) {
    if(data.data[j] > ymax || isnan(ymax))  // Stokes I
      ymax = data.data[j];
    if(data.data[j] < ymin || isnan(ymin))
      ymin = data.data[j];
    if(data.data[j+data.NrBins] > ymax)     // L
      ymax = data.data[j+data.NrBins];
    if(data.data[j+data.NrBins] < ymin)
      ymin = data.data[j+data.NrBins];
    if(data.data[j+2*data.NrBins] > ymax)   // V
      ymax = data.data[j+2*data.NrBins];
    if(data.data[j+2*data.NrBins] < ymin)
      ymin = data.data[j+2*data.NrBins];
    if(showtotpol) {
      if(data.data[j+5*data.NrBins] > ymax)    // Total polarization, but only check if you want to plot this
	ymax = data.data[j+5*data.NrBins];
      if(data.data[j+5*data.NrBins] < ymin)
	ymin = data.data[j+5*data.NrBins];
    }
  }

  //  printerror(verbose.debug, "XXXXXX %f %f", ymin, ymax);

  // Profile panel
  //  ppgsvp(0.1+pgplot->viewport.dxplot, 0.1+pgplot->viewport.dxplot+0.8*pgplot->viewport.xsize*overallplotscaling, 0.35+(1-overallplotscaling)*0.6+pgplot->viewport.dyplot, 0.35+(1-overallplotscaling)*0.6+pgplot->viewport.dyplot+0.55*pgplot->viewport.ysize*overallplotscaling);
  ppgsvp(vp_left, vp_right, vp_profile_bottom, vp_profile_top);

  frame.swin_x1 = longitude_left;
  frame.swin_x2 = longitude_right;
  if(xunit_type == 1) {
    frame.swin_x1 /= 360.0;
    frame.swin_x2 /= 360.0; 
  }
  if(Imin != Imax) {
    frame.swin_y1 = Imin;
    frame.swin_y2 = Imax;
    ppgswin(frame.swin_x1,frame.swin_x2,Imin,Imax);
  }else {
    frame.swin_y1 = ymin-(ymax-ymin)*0.05;
    frame.swin_y2 = ymax+(ymax-ymin)*0.05;
    ppgswin(frame.swin_x1,frame.swin_x2,ymin-(ymax-ymin)*0.05,ymax+(ymax-ymin)*0.05);
  }


  if(nopaswing == 1 && showEll == 0 && showPAdist == 0 && showELLdist == 0) {
    strcpy(pgplot->box.xlabel, xlabel);
    strcpy(pgplot->box.box_xopt, "bcnst");
  }else {
    pgplot->box.xlabel[0] = 0;
    strcpy(pgplot->box.box_xopt, "bcst");
  }

  //  ppgslw(boxlinewidth);
  if(noynumbers) {
    strcpy(pgplot->box.box_yopt, "bcts");
    //    ppgbox("bcst",xtick,nxsub,"bcts",0.0,0);
  }else {
    strcpy(pgplot->box.box_yopt, "bcnts");
    //    ppgbox("bcst",xtick,nxsub,"bcnts",0.0,0);
  }
  //  ppgsch(1);
  if(!noynumbers) {
    pgplot->box.drawlabels = 1;
    strcpy(pgplot->box.ylabel, ylabel);
    //    ppglab("", ylabel, "");
  }
  //  ppgsch(titlech);
  //  ppgslw(titlelw);
  //  ppgmtxt("T", 1.0, 0.5, 0.5, title);
  //  ppgsch(1);

  pgplot_drawbox(&(pgplot->box));
  pgplot->box.drawtitle = 0;   /* Don't draw title above PA plot */

  ppgslw(datalinewidth);


  // Draw the Stokes V line
  if(dashed) {
    ppgsls(4);
  }else
    ppgsci(3);
  domove = 1;
  for(j = 0; j < (data.NrBins); j++) {
    float xpos;
    xpos = get_pulse_longitude(data, 0, j, verbose);
    if(xunit_type == 1) {
      xpos /= 360.0;
    }
    if(isnan(xpos) || isnan(data.data[j+2*data.NrBins])) {
      domove = 1;
    }else {
      if(domove) {
	ppgmove(xpos, data.data[j+2*data.NrBins]);
	domove = 0;
      }else {
	ppgdraw(xpos, data.data[j+2*data.NrBins]);
	//	printerror(verbose.debug, "XXXXX %f %f", get_pulse_longitude(data, 0, j, verbose), data.data[j+2*data.NrBins]);
      }
    }
  }

  // Draw the L line
  if(dashed)
    ppgsls(2);
  else
    ppgsci(2);
  domove = 1;
  for(j = 0; j < (data.NrBins); j++) {
    float xpos;
    xpos = get_pulse_longitude(data, 0, j, verbose);
    if(xunit_type == 1) {
      xpos /= 360.0;
    }
    if(isnan(xpos) || isnan(data.data[j+1*data.NrBins])) {
      domove = 1;
    }else {
      if(domove) {
	ppgmove(xpos, data.data[j+1*data.NrBins]);
	domove = 0;
      }else {
	ppgdraw(xpos, data.data[j+1*data.NrBins]);
	//	printerror(verbose.debug, "XXXXX %f %f", get_pulse_longitude(data, 0, j, verbose), data.data[j+2*data.NrBins]);
      }
    }
  }

  // Draw the total polarization
  if(showtotpol) {
    if(dashed) {
      ppgsls(1);
    }
    if(!dashed) {
      ppgsci(4);
      ppgsls(2);
    }

    //    ppgmove(get_pulse_longitude(data, 0, 0, verbose), Ppulse[0]);
    //    for(j = 1; j < (data.NrBins); j++) {
    //      ppgdraw(get_pulse_longitude(data, 0, j, verbose), Ppulse[j]);
    //    }

    domove = 1;
    for(j = 0; j < (data.NrBins); j++) {
      float xpos;
      xpos = get_pulse_longitude(data, 0, j, verbose);
      if(xunit_type == 1) {
	xpos /= 360.0;
      }
      if(isnan(xpos) || isnan(data.data[j+5*data.NrBins])) {
	domove = 1;
      }else {
	if(domove) {
	  ppgmove(xpos, data.data[j+5*data.NrBins]);
	  domove = 0;
	}else {
	  ppgdraw(xpos, data.data[j+5*data.NrBins]);
	  //	printerror(verbose.debug, "XXXXX %f %f", get_pulse_longitude(data, 0, j, verbose), data.data[j+5*data.NrBins]);
	}
      }
    }
  }

  // Draw Stokes I
  ppgsci(1);
  ppgsls(1);
  domove = 1;
  for(j = 0; j < (data.NrBins); j++) {
    float xpos;
    xpos = get_pulse_longitude(data, 0, j, verbose);
    if(xunit_type == 1) {
      xpos /= 360.0;
    }
    //    printerror(verbose.debug, "XXXXX %f %f", get_pulse_longitude(data, 0, j, verbose), data.data[j]);
    if(isnan(xpos) || isnan(data.data[j])) {
      domove = 1;
    }else {
      if(domove) {
	ppgmove(xpos, data.data[j]);
	domove = 0;
      }else {
	ppgdraw(xpos, data.data[j]);
	//	printerror(verbose.debug, "XXXXX %f %f", get_pulse_longitude(data, 0, j, verbose), data.data[j]);
      }
    }
  }


  ppgsls(1);

  if(textoption != NULL && argv != NULL && argc != 0) {
    for(j = 0; j < argc; j++) {
      if(strcmp(argv[j], textoption) == 0) {
	i = sscanf(argv[j+1], "%f %f %f %d %d %d", &text_x, &text_y, &text_ch, &text_lw, &text_f, &text_ci);
	if(xunit_type == 1) {
	  text_x /= 360.0;
	}
	if(i != 6) {
	  fflush(stdout);
	  printerror(verbose.debug, "ERROR pgplotPAplot: Error parsing %s option", textoption);
	  memcpy(pgplot, pgplot_backup, sizeof(pgplot_options_definition));
	  free(pgplot_backup);
	  return 0;
	}
	ppgsch(text_ch);
	ppgscf(text_f);
	newtext = str_replace_header_params(data, argv[j+2], verbose);
	if(newtext == NULL) {
	  printerror(verbose.debug, "ERROR pgplotPAplot: Cannot substuture header parameter in -text option");
	  memcpy(pgplot, pgplot_backup, sizeof(pgplot_options_definition));
	  free(pgplot_backup);
	  return 0;
	}
	if(outline_txt > 0) {
	  ppgslw(text_lw+outline_lw);
	  ppgsci(outline_color);
	  ppgptxt(text_x, text_y, 0, 0, newtext);
	}
	ppgsci(text_ci);
	ppgslw(text_lw);
	ppgptxt(text_x, text_y, 0, 0, newtext);
	ppgscf(1);
	ppgsch(1);
	ppgslw(1);
	ppgsci(1);
	j += 2;
	free(newtext);
      }
    }
  }

  /* Note: same code is used a little bit below */
  if(herrorbaroption != NULL && argv != NULL && argc != 0) {
    for(j = 0; j < argc; j++) {
      float herr_x1, herr_x2, herr_x3, herr_y, herr_size, herr_lw;
      int herr_ci;
      if(strcmp(argv[j], herrorbaroption) == 0) {
	i = sscanf(argv[j+1], "%f %f %f %f %f %f %d", &herr_x1, &herr_x2, &herr_x3, &herr_y, &herr_size, &herr_lw, &herr_ci);
	if(xunit_type == 1) {
	  herr_x1 /= 360.0;
	  herr_x2 /= 360.0;
	  herr_x3 /= 360.0;
	}
	if(i != 7) {
	  fflush(stdout);
	  printerror(verbose.debug, "ERROR pgplotPAplot: Error parsing %s option", herrorbaroption);
	  memcpy(pgplot, pgplot_backup, sizeof(pgplot_options_definition));
	  free(pgplot_backup);
	  return 0;
	}
	ppgsch(herr_size);
	ppgsci(herr_ci);
	ppgslw(herr_lw);
	ppgpt1(herr_x2, herr_y*(frame.swin_y2-frame.swin_y1)+frame.swin_y1, 2);
	ppgerr1(1, herr_x2, herr_y*(frame.swin_y2-frame.swin_y1)+frame.swin_y1, herr_x3-herr_x2, herr_size);
	ppgerr1(3, herr_x2, herr_y*(frame.swin_y2-frame.swin_y1)+frame.swin_y1, herr_x2-herr_x1, herr_size);
	ppgsch(1);
	ppgslw(1);
	ppgsci(1);
	j += 1;
      }
    }
  }

  /* Note: same code is used a little bit below */
  if(verrorbaroption != NULL && argv != NULL && argc != 0) {
    for(j = 0; j < argc; j++) {
      float verr_y1, verr_y2, verr_y3, verr_x, verr_size, verr_lw;
      int verr_ci;
      if(strcmp(argv[j], verrorbaroption) == 0) {
	i = sscanf(argv[j+1], "%f %f %f %f %f %f %d", &verr_x, &verr_y1, &verr_y2, &verr_y3, &verr_size, &verr_lw, &verr_ci);
	if(xunit_type == 1) {
	  verr_x /= 360.0;
	}
	if(i != 7) {
	  fflush(stdout);
	  printerror(verbose.debug, "ERROR pgplotPAplot: Error parsing %s option", verrorbaroption);
	  memcpy(pgplot, pgplot_backup, sizeof(pgplot_options_definition));
	  free(pgplot_backup);
	  return 0;
	}
	ppgsch(verr_size);
	ppgsci(verr_ci);
	ppgslw(verr_lw);
	ppgpt1(verr_x, verr_y2*(frame.swin_y2-frame.swin_y1)+frame.swin_y1, 2);
	ppgerr1(2, verr_x, verr_y2*(frame.swin_y2-frame.swin_y1)+frame.swin_y1, verr_y3-verr_y2, verr_size);
	ppgerr1(4, verr_x, verr_y2*(frame.swin_y2-frame.swin_y1)+frame.swin_y1, verr_y2-verr_y1, verr_size);
	ppgsch(1);
	ppgslw(1);
	ppgsci(1);
	j += 1;
      }
    }
  }

  // START OF PA-SWING CODE
  if(nopaswing == 0) {
    ppgsls(1);
    ppgsci(1);
    ppgslw(1);
    //  ppgsvp(0.1+pgplot->viewport.dxplot,0.1+pgplot->viewport.dxplot+0.8*pgplot->viewport.xsize*overallplotscaling, 0.35+(1-overallplotscaling)*0.6+pgplot->viewport.dyplot-0.25*pgplot->viewport.ysize*ysize2*overallplotscaling, 0.35+(1-overallplotscaling)*0.6+pgplot->viewport.dyplot);
    ppgsvp(vp_left, vp_right, vp_paswing_bottom, vp_paswing_top);
    
    frame.swin_x1 = longitude_left;
    frame.swin_x2 = longitude_right;
    if(xunit_type == 1) {
      frame.swin_x1 /= 360.0;
      frame.swin_x2 /= 360.0; 
    }
    frame.swin_y1 = pa_bottom;
    frame.swin_y2 = pa_top;
    ppgswin(frame.swin_x1,frame.swin_x2,pa_bottom,pa_top);
    //
    ppgslw(pgplot->box.box_lw);
    //  if(noynumbers) {
    //    ppgbox("bcnst",xtick, nxsub,"bcst",0.0,0);
    //  }else {
    //    ppgbox("bcnst",xtick, nxsub,"bcnst",0.0,0);
    //  }
    pgplot->box.drawlabels = 1;
    if(showEll == 0 && showPAdist == 0 && showELLdist == 0) {
      strcpy(pgplot->box.xlabel, xlabel);
      strcpy(pgplot->box.box_xopt, "bcnst");
    }else {
      pgplot->box.xlabel[0] = 0;
      strcpy(pgplot->box.box_xopt, "bcst");
    }
    if(!noynumbers) {
      strcpy(pgplot->box.ylabel, ylabel_pa);
      //    ppglab(xlabel, ylabel_pa, "");
    }else {
      pgplot->box.ylabel[0] = 0;
      //    ppglab(xlabel, "", "");
    }
    pgplot_drawbox(&(pgplot->box));
    
    float x, xold;
    ppgslw(datalinewidth);
    if(overlayPAfine)         // Draw 100 points between each data point, allowing smoother curves and avoid problems when sweep is too steep.
      overlayPAfine = 100;
    else
      overlayPAfine = 1;
    if(overlayPA) {
      ppgsci(2);
      for(j = 0; j < (data.NrBins); j++) {
	for(i = 0; i < overlayPAfine; i++) {
	  x = get_pulse_longitude(data, 0, j, verbose);
	  if(j == data.NrBins-1) {
	    x += (get_pulse_longitude(data, 0, j, verbose)-get_pulse_longitude(data, 0, j-1, verbose))*i/(float)overlayPAfine;
	  }else {
	    x += (get_pulse_longitude(data, 0, j+1, verbose)-get_pulse_longitude(data, 0, j, verbose))*i/(float)overlayPAfine;
	  }
	  I = paswing(overlayalpha, overlaybeta, x, overlaypa0, overlayl0, nrJumps, jump_longitudes, jump_offsets, 0, 0);
	  if(xunit_type == 1) {
	    x /= 360.0;
	  }
	  
	  if(I > 180)
	    I -= 180;
	  if(I > 180)
	    I -= 180;
	  if(I < 0)
	    I += 180;
	  if(I < 0)
	    I += 180;
	  if(j != 0 || i != 0) {
	    if(fabs(I-Iold-180) < fabs(I-Iold))
	      Iold -= 180;
	    if(fabs(I-Iold+180) < fabs(I-Iold))
	      Iold += 180;
	    if(fabs(I-Iold) < 30) {
	      ppgmove(xold, Iold);
	      ppgdraw(x, I);
	      if(pa_bottom < -10) {
		ppgmove(xold, Iold-180);
		ppgdraw(x, I-180);
	      }
	    }
	  }
	  Iold = I;
	  xold = x;
	}
      }
      ppgsci(1);
    }

    ppgsci(1);
    for(j = 0; j < (data.NrBins); j++) {
      ok = 1;
      if(data.offpulse_rms != NULL) {
	if(data.data[j+data.NrBins] < sigma_limit*data.offpulse_rms[1])
	  ok = 0;
      }
      if(data.data[j+4*data.NrBins] < 0)
	ok = 0;
      if(ok && !isnan(data.data[j+3*data.NrBins]) && !isnan(get_pulse_longitude(data, 0, j, verbose))) {
	I = derotate_180(data.data[j+3*data.NrBins]+PAoffset) + 180.0;
	float xpos;
	xpos = get_pulse_longitude(data, 0, j, verbose);
	if(xunit_type == 1) {
	  xpos /= 360.0;
	}
	ppgpt1(xpos, I, -2);
	//      ymax = I+data.data[j+4*data.NrBins];
	//      ymin = I-data.data[j+4*data.NrBins];
	ppgerr1(6, xpos, I, data.data[j+4*data.NrBins], 1.0);
	
	I = derotate_180(data.data[j+3*data.NrBins]+PAoffset);
	ppgpt1(xpos, I, -2);
	//      ymax = I+data.data[j+4*data.NrBins];
	//      ymin = I-data.data[j+4*data.NrBins];
	ppgerr1(6, xpos, I, data.data[j+4*data.NrBins], 1.0);
	
	I = derotate_180(data.data[j+3*data.NrBins]+PAoffset) - 180.0;
	ppgpt1(xpos, I, -2);
	//      ymax = I+data.data[j+4*data.NrBins];
	//      ymin = I-data.data[j+4*data.NrBins];
	ppgerr1(6, xpos, I, data.data[j+4*data.NrBins], 1.0);
	
	I = derotate_180(data.data[j+3*data.NrBins]+PAoffset) - 360.0;
	ppgpt1(xpos, I, -2);
	//      ymax = I+data.data[j+4*data.NrBins];
	//      ymin = I-data.data[j+4*data.NrBins];
	ppgerr1(6, xpos, I, data.data[j+4*data.NrBins], 1.0);


	/*
	  ppgpt1(xpos, data.data[j+3*data.NrBins]+PAoffset, -2); 
	  ymax = data.data[j+3*data.NrBins]+data.data[j+4*data.NrBins]+PAoffset;
	  ymin = data.data[j+3*data.NrBins]-data.data[j+4*data.NrBins]+PAoffset;
	  ppgerr1(6, xpos, data.data[j+3*data.NrBins]+PAoffset, data.data[j+4*data.NrBins], 1.0);
	  I = derotate_deg(data.data[j+3*data.NrBins]+180+PAoffset);
	  if(I > 180)
	  I -= 360;
	  ppgpt1(xpos, I, -2);
	  ymax = I+data.data[j+4*data.NrBins];
	  ymin = I-data.data[j+4*data.NrBins];
	  ppgerr1(6, xpos, I, data.data[j+4*data.NrBins], 1.0);
	*/
      }
    }

    /* Note: same code is used a little bit below */
    if(herrorbaroption2 != NULL && argv != NULL && argc != 0) {
      for(j = 0; j < argc; j++) {
	float herr_x1, herr_x2, herr_x3, herr_y, herr_size, herr_lw;
	int herr_ci;
	if(strcmp(argv[j], herrorbaroption2) == 0) {
	  i = sscanf(argv[j+1], "%f %f %f %f %f %f %d", &herr_x1, &herr_x2, &herr_x3, &herr_y, &herr_size, &herr_lw, &herr_ci);
	  if(xunit_type == 1) {
	    herr_x1 /= 360.0;
	    herr_x2 /= 360.0;
	    herr_x3 /= 360.0;
	  }
	  if(i != 7) {
	    fflush(stdout);
	    printerror(verbose.debug, "ERROR pgplotPAplot: Error parsing %s option", herrorbaroption2);
	    memcpy(pgplot, pgplot_backup, sizeof(pgplot_options_definition));
	    free(pgplot_backup);
	    return 0;
	  }
	  ppgsch(herr_size);
	  ppgsci(herr_ci);
	  ppgslw(herr_lw);
	  ppgpt1(herr_x2, herr_y, 2);
	  ppgerr1(1, herr_x2, herr_y, herr_x3-herr_x2, herr_size);
	  ppgerr1(3, herr_x2, herr_y, herr_x2-herr_x1, herr_size);
	  ppgsch(1);
	  ppgslw(1);
	  ppgsci(1);
	  j += 1;
	}
      }
    }
    
    /* Note: same code is used a little bit below */
    if(verrorbaroption2 != NULL && argv != NULL && argc != 0) {
      for(j = 0; j < argc; j++) {
	float verr_y1, verr_y2, verr_y3, verr_x, verr_size, verr_lw;
	int verr_ci;
	if(strcmp(argv[j], verrorbaroption2) == 0) {
	  i = sscanf(argv[j+1], "%f %f %f %f %f %f %d", &verr_x, &verr_y1, &verr_y2, &verr_y3, &verr_size, &verr_lw, &verr_ci);
	  if(xunit_type == 1) {
	    verr_x /= 360.0;
	  }
	  if(i != 7) {
	    fflush(stdout);
	    printerror(verbose.debug, "ERROR pgplotPAplot: Error parsing %s option", verrorbaroption);
	    memcpy(pgplot, pgplot_backup, sizeof(pgplot_options_definition));
	    free(pgplot_backup);
	    return 0;
	  }
	  ppgsch(verr_size);
	  ppgsci(verr_ci);
	  ppgslw(verr_lw);
	  ppgpt1(verr_x, verr_y2, 2);
	  ppgerr1(2, verr_x, verr_y2, verr_y3-verr_y2, verr_size);
	  ppgerr1(4, verr_x, verr_y2, verr_y2-verr_y1, verr_size);
	  ppgsch(1);
	  ppgslw(1);
	  ppgsci(1);
	  j += 1;
	}
      }
    }
  }

  // START OF ELLIPTICITY STUFF
  if(showEll) {
    ppgsls(1);
    ppgsci(1);
    ppgslw(1);
    //  ppgsvp(0.1+pgplot->viewport.dxplot, 0.1+pgplot->viewport.dxplot+0.8*pgplot->viewport.xsize*overallplotscaling, 0.35+(1-overallplotscaling)*0.6+pgplot->viewport.dyplot, 0.35+(1-overallplotscaling)*0.6+pgplot->viewport.dyplot+0.55*pgplot->viewport.ysize*overallplotscaling);
    //   ppgsvp(0.1+pgplot->viewport.dxplot,0.1+pgplot->viewport.dxplot+0.8*pgplot->viewport.xsize*overallplotscaling, 0.35+(1-overallplotscaling)*0.6+pgplot->viewport.dyplot-0.25*pgplot->viewport.ysize*ysize2*overallplotscaling, 0.35+(1-overallplotscaling)*0.6+pgplot->viewport.dyplot);
    //    ppgsvp(0.1+pgplot->viewport.dxplot, 0.1+pgplot->viewport.dxplot+0.8*pgplot->viewport.xsize*overallplotscaling, 0.35+(1-overallplotscaling)*(0.6 + 0.25*pgplot->viewport.ysize*ysize2)+pgplot->viewport.dyplot-0.25*pgplot->viewport.ysize*ysize2-0.25*pgplot->viewport.ysize*ysize2*overallplotscaling, 0.35+(1-overallplotscaling)*(0.6 + 0.25*pgplot->viewport.ysize*ysize2)+pgplot->viewport.dyplot-0.25*pgplot->viewport.ysize*ysize2);
    ppgsvp(vp_left, vp_right, vp_ellswing_bottom, vp_ellswing_top);

    frame.swin_x1 = longitude_left;
    frame.swin_x2 = longitude_right;
    if(xunit_type == 1) {
      frame.swin_x1 /= 360.0;
      frame.swin_x2 /= 360.0;
    }
    frame.swin_y1 = -45;
    frame.swin_y2 = +45;
    ppgswin(frame.swin_x1,frame.swin_x2,-45,45);
    ppgslw(pgplot->box.box_lw);
    //  if(noynumbers) {
    //    ppgbox("bcnst",xtick, nxsub,"bcst",0.0,0);
    //  }else {
    //    ppgbox("bcnst",xtick, nxsub,"bcnst",0.0,0);
    //  }
    pgplot->box.drawlabels = 1;
    if(!noynumbers) {
      strcpy(pgplot->box.ylabel, ylabel_ell);
    }else {
      pgplot->box.ylabel[0] = 0;
    }
    if(showPAdist == 0 && showELLdist == 0) {
      strcpy(pgplot->box.xlabel, xlabel);
      strcpy(pgplot->box.box_xopt, "bcnst");
    }else {
      pgplot->box.xlabel[0] = 0;
      strcpy(pgplot->box.box_xopt, "bcst");
    }
    pgplot_drawbox(&(pgplot->box));

    //    float x, xold;
    ppgslw(datalinewidth);

    for(j = 0; j < (data.NrBins); j++) {
      ok = 1;
      if(data.offpulse_rms != NULL) {
	if(data.data[j+5*data.NrBins] < sigma_limit*data.offpulse_rms[5])
	  ok = 0;
      }
      if(data.data[j+7*data.NrBins] < 0)
	ok = 0;
      if(ok && !isnan(data.data[j+6*data.NrBins]) && !isnan(get_pulse_longitude(data, 0, j, verbose))) {
	//	ppgsci(3);
	
	//	I = derotate_90(data.data[j+6*data.NrBins]) + 180.0;
	//	ppgpt1(get_pulse_longitude(data, 0, j, verbose), I, -2);
	//	ppgerr1(6, get_pulse_longitude(data, 0, j, verbose), I, data.data[j+7*data.NrBins], 1.0);
	
	I = derotate_90(data.data[j+6*data.NrBins]) + 90.0;
	float xpos;
	xpos = get_pulse_longitude(data, 0, j, verbose);
	if(xunit_type == 1) {
	  xpos /= 360.0;
	}
	ppgpt1(xpos, I, -2);
	ppgerr1(6, xpos, I, data.data[j+7*data.NrBins], 1.0);
	
	I = derotate_90(data.data[j+6*data.NrBins]);
	ppgpt1(xpos, I, -2);
	ppgerr1(6, xpos, I, data.data[j+7*data.NrBins], 1.0);
	
	I = derotate_90(data.data[j+6*data.NrBins]) - 90.0;
	ppgpt1(xpos, I, -2);
	ppgerr1(6, xpos, I, data.data[j+7*data.NrBins], 1.0);
	
	I = derotate_90(data.data[j+6*data.NrBins]) - 180.0;
	ppgpt1(xpos, I, -2);
	ppgerr1(6, xpos, I, data.data[j+7*data.NrBins], 1.0);

	//	I = derotate_90(data.data[j+6*data.NrBins]) - 270.0;
	//	ppgpt1(xpos, I, -2);
	//	ppgerr1(6, xpos, I, data.data[j+7*data.NrBins], 1.0);
      }

    }
  }   // End of if(showEll)


  int cmaptype;
  cmaptype = pgplot_device_type(NULL, verbose);
  //      printf("cmaptype=%d\n", cmaptype);
  if(cmaptype <= 2)
    cmaptype = PPGPLOT_GRAYSCALE;
  else
    cmaptype = PPGPLOT_INVERTED_GRAYSCALE;

  if(showPAdist) {
    ppgsls(1);
    ppgsci(1);
    ppgslw(1);
    ppgsvp(vp_left, vp_right, vp_padist_bottom, vp_padist_top);

    //    frame.swin_x1 = longitude_left;
    //    frame.swin_x2 = longitude_right;
    //    frame.swin_y1 = -45;
    //    frame.swin_y2 = +45;
    //    ppgswin(longitude_left,longitude_right,-45,45);
    ppgslw(pgplot->box.box_lw);
    pgplot->box.drawlabels = 1;
    if(!noynumbers) {
      strcpy(pgplot->box.ylabel, ylabel_pa);
    }else {
      pgplot->box.ylabel[0] = 0;
    }
    if(showELLdist == 0) {
      strcpy(pgplot->box.xlabel, xlabel);
      strcpy(pgplot->box.box_xopt, "bcnst");
    }else {
      pgplot->box.xlabel[0] = 0;
      strcpy(pgplot->box.box_xopt, "bcst");
    }

    //    float x, xold;
    ppgslw(datalinewidth);
    pgplot->viewport.dontopen = 1;
    pgplot->viewport.dontclose = 1; 
    float xpos1, xpos2, xpos_start, xpos_end;
    xpos1 = get_pulse_longitude(data, 0, 0, verbose);
    xpos2 = get_pulse_longitude(data, 0, padist->NrBins-1, verbose);
    xpos_start = longitude_left;
    xpos_end = longitude_right;
    if(xunit_type) {
      xpos1 /= 360.0;
      xpos2 /= 360.0;
      xpos_start /= 360.0;
      xpos_end /= 360.0;
    }
    float ymin_data, ymax_data;
    ymin_data = -90;        // The nominal range of the data
    ymax_data = 90;
    //    printf("XXXX Range of input data %f %f\n", ymin_data, ymax_data);
    long wraps;
    wraps = floor((padist_pamin-ymin_data)/180.0);
    //    printf("XXXX Add %ld 180 deg wraps to cover requested lower bound\n", wraps);
    ymin_data += 180.0*wraps;
    ymax_data += 180.0*wraps;
    //    printf("XXXX New assigned range of data %f %f\n", ymin_data, ymax_data);
    wraps = ceil((padist_pamax-ymax_data)/180.0)+1;
    if(wraps < 1)
      wraps = 1;
    //    printf("XXXX Need to plot data range %ld times\n", wraps);
    float *data_ptr;
    if(wraps == 1) {
      data_ptr = padist->data;
    }else if(wraps > 10) {  // This should already be caught at the start of the function
      printerror(verbose.debug, "ERROR pgplotPAplot: Cannot plot pa-distribution more than 10 times above each other. Something appears to be wrong with the requested PA range for the PA distribution panel.\n");
      return 0;      
    }else {
      long blocksize;
      blocksize = padist->NrBins*padist->NrSubints;
      data_ptr = (float *)malloc(wraps*blocksize*sizeof(float));
      if(data_ptr == NULL) {
	printerror(verbose.debug, "ERROR pgplotPAplot: Memory allocation error\n");
	return 0;      
      }
      for(i = 0; i < wraps; i++) {
	memcpy(data_ptr+i*blocksize, padist->data, blocksize*sizeof(float));
      }
      ymax_data += 180.0*(wraps-1);
    }
    //    printf("XXXX New assigned range of data %f %f\n", ymin_data, ymax_data);

    ymin_data += 0.5*(180.0/(float)padist->NrSubints);
    ymax_data -= 0.5*(180.0/(float)padist->NrSubints);
    //    printf("XXXX New assigned range of data %f %f (correcting to bin centres)\n", ymin_data, ymax_data);
    if(pgplotMap(pgplot, data_ptr, padist->NrBins, wraps*padist->NrSubints, xpos1, xpos2, xpos_start, xpos_end, ymin_data, ymax_data, padist_pamin, padist_pamax, cmaptype, 0, 0, 0, NULL, 0, 0, padist_saturize, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, showwedge, 1, 0, verbose) == 0) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR pgplotPAplot: Plotting of the PA distribution failed.");
      memcpy(pgplot, pgplot_backup, sizeof(pgplot_options_definition));
      free(pgplot_backup);
      return 0;
    }
    if(wraps > 1) {
      free(data_ptr);
    }
    pgplot_drawbox(&(pgplot->box));

  }   // End of if(showPAdist)

  if(showELLdist) {
    ppgsls(1);
    ppgsci(1);
    ppgslw(1);
    ppgsvp(vp_left, vp_right, vp_elldist_bottom, vp_elldist_top);

    //    frame.swin_x1 = longitude_left;
    //    frame.swin_x2 = longitude_right;
    //    frame.swin_y1 = -45;
    //    frame.swin_y2 = +45;
    //    ppgswin(longitude_left,longitude_right,-45,45);
    ppgslw(pgplot->box.box_lw);
    pgplot->box.drawlabels = 1;
    strcpy(pgplot->box.xlabel, xlabel);
    if(!noynumbers) {
      strcpy(pgplot->box.ylabel, ylabel_ell);
    }else {
      pgplot->box.ylabel[0] = 0;
    }
    strcpy(pgplot->box.box_xopt, "bcnst");

    //    float x, xold;
    ppgslw(datalinewidth);
    pgplot->viewport.dontopen = 1;
    pgplot->viewport.dontclose = 1; 
    float xpos1, xpos2, xpos_start, xpos_end;
    xpos1 = get_pulse_longitude(data, 0, 0, verbose);
    xpos2 = get_pulse_longitude(data, 0, elldist->NrBins-1, verbose);
    xpos_start = longitude_left;
    xpos_end = longitude_right;
    if(xunit_type) {
      xpos1 /= 360.0;
      xpos2 /= 360.0;
      xpos_start /= 360.0;
      xpos_end /= 360.0;
    }
    if(pgplotMap(pgplot, elldist->data, elldist->NrBins, elldist->NrSubints, xpos1, xpos2, xpos_start, xpos_end, -45+0.5*(90.0/(float)elldist->NrSubints), 45-0.5*(90.0/(float)elldist->NrSubints), -45, 45, cmaptype, 0, 0, 0, NULL, 0, 0, elldist_saturize, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, showwedge, 1, 0, verbose) == 0) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR pgplotPAplot: Plotting of the ellipticity angle distribution failed.");
      memcpy(pgplot, pgplot_backup, sizeof(pgplot_options_definition));
      free(pgplot_backup);
      return 0;
    }
    pgplot_drawbox(&(pgplot->box));

  }   // End of if(showElldist)

  if(pgplot->viewport.dontclose == 0)
    ppgclos();
  memcpy(pgplot, pgplot_backup, sizeof(pgplot_options_definition));
  free(pgplot_backup);
  return 1;
}


//START REGION DEVELOP
//START REGION RELEASE

/*
  Plot 'data' with nrx points. Properties of the device is set with
  viewport and properties of the box and labels with pgplotbox. If
  datax is set, those x values will be used. If datax is set to NULL,
  the data should be uniformely sampled between xmin and xmax, so data
  is an one dimensional array. If sigma != NULL, it will be used to
  plot errorbars. If ymin_show == ymax_show the minumum and maximum
  are obtained from the data (if forceMinZero is set than the yrange
  is set from 0 to the maximum in the data). If ymin_show !=
  ymax_show, these will used for the vertical range. Set xmin_show and
  xmax_show to zoom (if they are equal the full range is used
  instead). If dontsetranges is set, the ranges are expected to be
  already by previous pgplot calls. If hist is set, a histogram is
  produced rather than connected points (if forceMinZero is set, it
  means the outer edges of the histogram will go to zero.). If noline
  is set there is no line drawn in between the points. The width of
  the curve can be set with linewidth. If pointtype is set the points
  are drawn with this pgplot symbol and color sets the color. boxcolor
  sets the color of the axis (default = 1). You can specify regions by
  using the regions struct, or you can set it to NULL if not
  used. This will increase the color index from color (offpulse) by 1
  for each subsequent selected onpulse region. Alternatively, if
  onpulsecolor is non-negative, the color is onpulsecolor inside the
  selected region(s) and color outside the regions.

  Return 1: Succesful
         0: Error
*/

int pgplotGraph1(pgplot_options_definition *pgplot, float *data, float *datax, float *sigma, int nrx, float xmin, float xmax, int dontsetranges, float xmin_show, float xmax_show, float ymin_show, float ymax_show, int forceMinZero, int hist, int noline, int linewidth, int pointtype, int color, int boxcolor, pulselongitude_regions_definition *regions, int onpulsecolor, verbose_definition verbose)
{
  float min, max;
  float x, xnext, y;
  int binnr, regnr, loopnr, deviceID, notset, color_prev, color_cur, curpencolor;
  pgplot_frame_def_internal frame;

  if(hist && datax != NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR pgplotGraph1: Cannot plot histogram when x-values are provided");
    return 0;    
  }
  if(hist && nrx == 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR pgplotGraph1: Cannot plot histogram for one point");
    return 0;    
  }
  if(hist == 0 && regions != NULL) {
    fflush(stdout);
    printwarning(verbose.debug, "WARNING pgplotGraph1: When showing selected regions, histogram mode is more precise");
  }
  if(regions != NULL) {
    for(regnr = 0; regnr < regions->nrRegions; regnr++) {
      if(regions->bins_defined[regnr] == 0) {
	fflush(stdout);
	printerror(verbose.debug, "ERROR pgplotGraph1: onpulse region is not defined in bins");
	return 0;
      }
    }
  }

  // Shouldn't be necessary, but gets rid of warning when compiling with optimisation
  xnext = 0;

  //  fflush(stdout); fprintf(stderr, "XXXXXX open device\n");
  if(pgplot_opendevice(&(pgplot->viewport), &deviceID, verbose) == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR pgplotGraph1: Cannot open plot device");
    return 0;
  }

  //  ppgslw(1);
  //  ppgsvp(0.1+viewport->dxplot, 0.1+viewport->dxplot+0.8*viewport->xsize, 0.1+viewport->dyplot, 0.1+viewport->dyplot+0.8*viewport->ysize);

  notset = 1;
  min = max = 0;
  if(datax != NULL) {
    xmin = xmax = datax[0];
  }
  for(binnr = 0; binnr < nrx; binnr++) {
    if(datax == NULL) {
      if(nrx == 1) {
	x = xmin;
      }else {
	x = binnr*(xmax-xmin)/(float)(nrx-1) + xmin;
      }
    }else {
      x = datax[binnr];
    }
    if(xmin_show == xmax_show || (xmin_show < xmax_show && x >= xmin_show && x <= xmax_show) || (xmin_show > xmax_show && x <= xmin_show && x >= xmax_show)) {
      if(notset) {
	min = max = data[binnr];
	notset = 0;
      }
      if(data[binnr] > max)
	max = data[binnr];
      if(data[binnr] < min)
	min = data[binnr];
    }
    if(datax != NULL) {
      if(datax[binnr] > xmax)
	xmax = datax[binnr];
      if(datax[binnr] < xmin)
	xmin = datax[binnr];
    }
  }
  if(forceMinZero)
    min = 0;

  if(min == max) {
    min -= 1.0;
    max += 1.0;
  }

  if(ymin_show != ymax_show) {
    min = ymin_show;
    max = ymax_show;
  }else {
    if(forceMinZero == 0 || hist == 0)
      min -= (max-min)*0.01;
    max += (max-min)*0.01;
  }

  if(xmin_show == xmax_show) {
    //    printf("XXXXX  Setting %f %f\n", xmin, xmax);
    xmin_show = xmin;
    xmax_show = xmax;
  }

  if(hist) {
    xmin_show -= 0.5*(xmax-xmin)/(float)(nrx-1);
    xmax_show += 0.5*(xmax-xmin)/(float)(nrx-1);
  }

  //  fflush(stdout); fprintf(stderr, "XXXXXX clear frame\n");
  clear_pgplot_frame(&frame);
  frame.svp = 1;
  frame.svp_x1 = 0.1+pgplot->viewport.dxplot;
  frame.svp_x2 = 0.1+pgplot->viewport.dxplot+0.8*pgplot->viewport.xsize;
  frame.svp_y1 = 0.1+pgplot->viewport.dyplot;
  frame.svp_y2 = 0.1+pgplot->viewport.dyplot+0.8*pgplot->viewport.ysize;
  if(dontsetranges == 0)
    frame.swin = 1;
  else
    frame.swin = 0;
  frame.swin_x1 = xmin_show;
  frame.swin_x2 = xmax_show;
  frame.swin_y1 = min;
  frame.swin_y2 = max;
  //  ppgswin(xmin_show, xmax_show, min, max);
  //  fflush(stdout); fprintf(stderr, "XXXXXX make frame\n");
  pgplot_makeframe(&frame);

 //  ppgbox("bcnsti",0.0,0,"bcntsi",0.0,0);
  //  ppglab(xlabel, ylabel, title);
  ppgsci(boxcolor);
  curpencolor = boxcolor;
  //  fflush(stdout); fprintf(stderr, "XXXXXX draw box\n");
  pgplot_drawbox(&(pgplot->box));


  color_prev = color;

  //  fflush(stdout); fprintf(stderr, "XXXXXX start buffering\n");
  ppgbbuf();
  ppgslw(linewidth);
  for(loopnr = 0; loopnr < 2; loopnr++) {                    /* If k=0, draw line, if k=1, draw points. */
    for(binnr = 0; binnr < nrx; binnr++) {
      if(datax == NULL) {
	if(nrx == 1) {
	  x = xmin;
	  xnext = xmin;
	}else {
	  x     = binnr*(xmax-xmin)/(float)(nrx-1) + xmin;
	  xnext = (binnr+1)*(xmax-xmin)/(float)(nrx-1) + xmin;
	}
	if(hist) {
	  x -= 0.5*(xmax-xmin)/(float)(nrx-1);
	  xnext -= 0.5*(xmax-xmin)/(float)(nrx-1);
	}
      }else {
	x = datax[binnr];
      }
      y = data[binnr];
      color_cur = color;
      if(regions != NULL) {
	for(regnr = 0; regnr < regions->nrRegions; regnr++) {
	  //	  printf("XXXXXX regions: %d %d\n", regions->left_bin[regnr], regions->right_bin[regnr]);
	  if(hist == 0) {
	    if(binnr > regions->left_bin[regnr] && binnr <= regions->right_bin[regnr]) {
	      if(onpulsecolor < 0) {
		color_cur = color+1+regnr;
	      }else {
		color_cur = onpulsecolor;
	      }
	    }
	  }else {
	    if(binnr >= regions->left_bin[regnr] && binnr <= regions->right_bin[regnr]) {
	      if(onpulsecolor < 0) {
		color_cur = color+1+regnr;
	      }else {
		color_cur = onpulsecolor;
	      }
	    }
	  }
	}
      }
      if(loopnr == 0) {
	if(noline == 0) {
	  if(hist == 0) {
	    if(binnr == 0) {
	      ppgmove(x, y);
	    }else {
	      if(color_cur != curpencolor) {
		ppgsci(color_cur);
		curpencolor = color_cur;
	      }
	      ppgdraw(x, y);
	    }
	  }else { // Drawing histogram
	    if(binnr == 0) {
	      if(forceMinZero) { // If set, histogram starts and ends at y=0
		ppgmove(x, 0);
		if(color_cur != curpencolor) {
		  ppgsci(color_cur);
		  curpencolor = color_cur;
		}
		ppgdraw(x, y);
	      }else {
		ppgmove(x, y);
	      }
	    }else {
	      if(color_cur == color || color_prev == color) {   // Transitions to/from unselected region: make default color
		if(color != curpencolor) {
		  ppgsci(color);
		  curpencolor = color;
		}
	      }else {
		if(color_cur == color_prev) {
		  if(color_cur != curpencolor) {
		    ppgsci(color_cur);
		    curpencolor = color_cur;
		  }
		}else {
		  if(color_prev != curpencolor) {   // Make left side bin color of prev bin
		    ppgsci(color_prev);
		    curpencolor = color_prev;
		  }
		}
	      }
	      ppgdraw(x, y);
	    }
	    ppgsci(color_cur);
	    ppgdraw(xnext, y);
	    if(binnr == nrx - 1) {
	      if(forceMinZero) { // If set, histogram starts and ends at y=0
		ppgdraw(xnext, 0);
	      }
	    }
	  }
	}
      }else {
	ppgsci(color);
	if(pointtype) {
	  ppgpt1(x, y, pointtype);
	}
	if(sigma != NULL) {
	  ppgerr1(6, x, y, sigma[binnr], 1.0);
	}
      }
      color_prev = color_cur;
    }
  }
  //  fflush(stdout); fprintf(stderr, "XXXXXX end buffering\n");
  ppgebuf();
  //  fflush(stdout); fprintf(stderr, "XXXXXX end buffering done\n");

  ppgsci(1);
  if(!pgplot->viewport.dontclose)
    ppgclos();
  //  fflush(stdout); fprintf(stderr, "XXXXXX close device\n");
  return 1;
}

//START REGION DEVELOP
//START REGION RELEASE

/* Opens a pgplot device (deviceID is set) and sets the resolution if
   provided. The viewport is also cleared, unless the noclear flag is
   set. Returns 0 if not successful. */
int pgplot_opendevice(pgplot_viewport_definition *viewport, int *deviceID, verbose_definition verbose)
{
  if(viewport->dontopen == 0) {
    *deviceID = ppgopen(viewport->plotDevice);
    if(*deviceID <= 0) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR pgplot_opendevice: Cannot open plot device called '%s'", viewport->plotDevice);
      return 0;
    }
    pgplot_setWindowsize(viewport->windowwidth, viewport->windowheight, viewport->aspectratio);
    ppgask(0);
  }
  /* Needs to be outside dontopen flag, as sometimes the same device is used for multiple plots. */
  if(viewport->noclear == 0)
    ppgpage();
  return 1;
}

//START REGION DEVELOP
//START REGION RELEASE

void pgplot_clear_options(pgplot_options_definition *pgplot) {
  pgplot_clear_viewport_def(&(pgplot->viewport));
  clear_pgplot_box(&(pgplot->box));
}

/* Sets viewport parameters to default values. */
void pgplot_clear_viewport_def(pgplot_viewport_definition *viewport)
{
  viewport->windowwidth = 0;
  viewport->windowheight = 0;
  viewport->aspectratio = -1;
  viewport->dxplot = 0;
  viewport->xsize = 1; 
  viewport->dyplot = 0;
  viewport->ysize = 1; 
  viewport->noclear = 0;
  viewport->dontopen = 0;
  viewport->dontclose = 0;
  sprintf(viewport->plotDevice, "?");
}

//START REGION DEVELOP
//START REGION RELEASE

/*
  Returns 0 if there was an error.

  The color maps are:
    PPGPLOT_GRAYSCALE
    PPGPLOT_INVERTED_GRAYSCALE
    etc.
*/

int pgplot_set_maptype(int maptype, verbose_definition verbose)
{
  float contrast, brightness;
  contrast = 1.0;
  brightness = 0.5;
  if(maptype == PPGPLOT_GRAYSCALE) {
    float datal[] = {0.0, 1.0};
    float datar[] = {1.0, 0.0};
    float datag[] = {1.0, 0.0};
    float datab[] = {1.0, 0.0};
    ppgctab(datal, datar, datag, datab, 2, contrast, brightness);
  }else if(maptype == PPGPLOT_INVERTED_GRAYSCALE) {
    float datal[] = {0.0, 1.};
    float datar[] = {0.0, 1.};
    float datag[] = {0.0, 1.};
    float datab[] = {0.0, 1.};
    ppgctab(datal, datar, datag, datab, 2, contrast, brightness);
  }else if(maptype == PPGPLOT_RED) {
    float datal[] = {0.0, 1.0};
    float datar[] = {1.0, 0.0};
    float datag[] = {0.0, 0.0};
    float datab[] = {0.0, 0.0};
    ppgctab(datal, datar, datag, datab, 2, contrast, brightness);
  }else if(maptype == PPGPLOT_INVERTED_RED) {
    float datal[] = {0.0, 1.0};
    float datar[] = {0.0, 1.0};
    float datag[] = {0.0, 0.0};
    float datab[] = {0.0, 0.0};
    ppgctab(datal, datar, datag, datab, 2, contrast, brightness);
  }else if(maptype == PPGPLOT_GREEN) {
    float datal[] = {0.0, 1.0};
    float datar[] = {0.0, 0.0};
    float datag[] = {1.0, 0.0};
    float datab[] = {0.0, 0.0};
    ppgctab(datal, datar, datag, datab, 2, contrast, brightness);
  }else if(maptype == PPGPLOT_INVERTED_GREEN) {
    float datal[] = {0.0, 1.0};
    float datar[] = {0.0, 0.0};
    float datag[] = {0.0, 1.0};
    float datab[] = {0.0, 0.0};
    ppgctab(datal, datar, datag, datab, 2, contrast, brightness);
  }else if(maptype == PPGPLOT_BLUE) {
    float datal[] = {0.0, 1.0};
    float datar[] = {0.0, 0.0};
    float datag[] = {0.0, 0.0};
    float datab[] = {1.0, 0.0};
    ppgctab(datal, datar, datag, datab, 2, contrast, brightness);
  }else if(maptype == PPGPLOT_INVERTED_BLUE) {
    float datal[] = {0.0, 1.0};
    float datar[] = {0.0, 0.0};
    float datag[] = {0.0, 0.0};
    float datab[] = {0.0, 1.0};
    ppgctab(datal, datar, datag, datab, 2, contrast, brightness);
  }else if(maptype == PPGPLOT_CYAN) {
    float datal[] = {0.0, 1.0};
    float datar[] = {0.0, 0.0};
    float datag[] = {1.0, 0.0};
    float datab[] = {1.0, 0.0};
    ppgctab(datal, datar, datag, datab, 2, contrast, brightness);
  }else if(maptype == PPGPLOT_INVERTED_CYAN) {
    float datal[] = {0.0, 1.0};
    float datar[] = {0.0, 0.0};
    float datag[] = {0.0, 1.0};
    float datab[] = {0.0, 1.0};
    ppgctab(datal, datar, datag, datab, 2, contrast, brightness);
  }else if(maptype == PPGPLOT_HEAT) {
    float datal[] = {0.0, 0.4, 0.6, 0.8, 1.0};
    float datar[] = {1.0, 1.0, 1.0, 0.5, 0.0};
    float datag[] = {1.0, 1.0, 0.5, 0.0, 0.0};
    float datab[] = {1.0, 0.3, 0.0, 0.0, 0.0};
    ppgctab(datal, datar, datag, datab, 5, contrast, brightness);
  }else if(maptype == PPGPLOT_INVERTED_HEAT) {
    float datal[] = {0.0, 0.2, 0.4, 0.6, 1.0};
    float datar[] = {0.0, 0.5, 1.0, 1.0, 1.0};
    float datag[] = {0.0, 0.0, 0.5, 1.0, 1.0};
    float datab[] = {0.0, 0.0, 0.0, 0.3, 1.0};
    ppgctab(datal, datar, datag, datab, 5, contrast, brightness);
  }else if(maptype == PPGPLOT_COLD) {
    float datal[] = {0.0, 0.4, 0.6, 0.8, 1.0};
    float datar[] = {1.0, 0.3, 0.0, 0.0, 0.0};
    float datag[] = {1.0, 1.0, 0.5, 0.0, 0.0};
    float datab[] = {1.0, 1.0, 1.0, 0.5, 0.0};
    ppgctab(datal, datar, datag, datab, 5, contrast, brightness);
  }else if(maptype == PPGPLOT_INVERTED_COLD) {
    float datal[] = {0.0, 0.2, 0.4, 0.6, 1.0};
    float datar[] = {0.0, 0.0, 0.0, 0.3, 1.0};
    float datag[] = {0.0, 0.0, 0.5, 1.0, 1.0};
    float datab[] = {0.0, 0.5, 1.0, 1.0, 1.0};
    ppgctab(datal, datar, datag, datab, 5, contrast, brightness);
//START REGION DEVELOP
  }else if(maptype == PPGPLOT_PLASMA) {
    float datal[] = {0.0, 0.4, 0.6, 0.8, 1.0};
    float datar[] = {1.0, 1.0, 0.5, 0.0, 0.0};
    float datag[] = {1.0, 0.3, 0.0, 0.0, 0.0};
    float datab[] = {1.0, 1.0, 1.0, 0.5, 0.0};
    ppgctab(datal, datar, datag, datab, 5, contrast, brightness);
  }else if(maptype == PPGPLOT_INVERTED_PLASMA) {
    float datal[] = {0.0, 0.2, 0.4, 0.6, 1.0};
    float datar[] = {0.0, 0.0, 0.5, 1.0, 1.0};
    float datag[] = {0.0, 0.0, 0.0, 0.3, 1.0};
    float datab[] = {0.0, 0.5, 1.0, 1.0, 1.0};
    ppgctab(datal, datar, datag, datab, 5, contrast, brightness);
  }else if(maptype == PPGPLOT_FOREST) {
    float datal[] = {0.0, 0.1, 0.3, 0.6, 1.0};
    float datar[] = {1.0, 0.3, 0.0, 0.0, 0.0};
    float datag[] = {1.0, 1.0, 1.0, 0.5, 0.0};
    float datab[] = {1.0, 1.0, 0.5, 0.0, 0.0};
    ppgctab(datal, datar, datag, datab, 5, contrast, brightness);
  }else if(maptype == PPGPLOT_INVERTED_FOREST) {
    float datal[] = {0.0, 0.4, 0.7, 0.9, 1.0};
    float datar[] = {0.0, 0.0, 0.0, 0.3, 1.0};
    float datag[] = {0.0, 0.5, 1.0, 1.0, 1.0};
    float datab[] = {0.0, 0.0, 0.5, 1.0, 1.0};
    ppgctab(datal, datar, datag, datab, 5, contrast, brightness);
  }else if(maptype == PPGPLOT_ALIEN_GLOW) {  // NOTE THAT THIS IS ALMOST THE SAME AS THE HEAT2 MAP
    float datal[] = {-0.7, 0.0, 0.17, 0.33, 0.50, 0.67, 0.83, 1.0, 1.5};
    float datar[] = { 1.0, 1.0, 1.0,  1.0,  0.6,  0.0,  0.0,  0.0, 0.0};
    float datag[] = { 1.0, 0.0, 0.6,  1.0,  1.0,  1.0,  0.0,  0.0, 0.0};
    float datab[] = { 1.0, 0.0, 0.0,  0.0,  0.3,  1.0,  0.8,  0.3, 0.0};
    ppgctab(datal, datar, datag, datab, 9, contrast, brightness);
  }else if(maptype == PPGPLOT_INVERTED_ALIEN_GLOW) {
    float datal[] = {-0.5, 0.0, 0.17, 0.33, 0.50, 0.67, 0.83, 1.0, 1.7};
    float datar[] = { 0.0, 0.0, 0.0,  0.0,  0.6,  1.0,  1.0,  1.0, 1.0};
    float datag[] = { 0.0, 0.0, 0.0,  1.0,  1.0,  1.0,  0.6,  0.0, 1.0};
    float datab[] = { 0.0, 0.3, 0.8,  1.0,  0.3,  0.0,  0.0,  0.0, 1.0};
    ppgctab(datal, datar, datag, datab, 9, contrast, brightness);
//START REGION RELEASE
  }else if(maptype == PPGPLOT_HEAT2) {
    float datal[] = {-0.5, 0.0, 0.17, 0.33, 0.50, 0.67, 0.83, 1.0, 1.7};
    float datar[] = { 1.0, 1.0, 1.0,  1.0,  0.6,  0.0,  0.0,  0.0, 0.0};
    float datag[] = { 1.0, 0.0, 0.6,  1.0,  1.0,  1.0,  0.0,  0.0, 0.0};
    float datab[] = { 1.0, 0.0, 0.0,  0.0,  0.3,  1.0,  0.8,  0.3, 0.0};
    ppgctab(datal, datar, datag, datab, 9, contrast, brightness);
  }else if(maptype == PPGPLOT_INVERTED_HEAT2) {
    float datal[] = {-0.5, 0.0, 0.17, 0.33, 0.50, 0.67, 0.83, 1.0, 1.7};
    float datar[] = { 0.0, 0.0,  0.0,  0.0,  0.6,  1.0,  1.0, 1.0, 1.0};
    float datag[] = { 0.0, 0.0,  0.0,  1.0,  1.0,  1.0,  0.6, 0.0, 1.0};
    float datab[] = { 0.0, 0.3,  0.8,  1.0,  0.3,  0.0,  0.0, 0.0, 1.0};
    ppgctab(datal, datar, datag, datab, 9, contrast, brightness);
  }else if(maptype == PPGPLOT_HEAT3) {
    float datal[] = {0.0, 0.2, 0.35, 0.5, 0.8, 1.0};
    float datar[] = {1.0, 1.0, 1.0, 1.0, 0.0, 0.0};
    float datag[] = {1.0, 1.0, 0.0, 0.0, 0.0, 0.0};
    float datab[] = {1.0, 0.0, 0.0, 1.0, 1.0, 0.0};
    ppgctab(datal, datar, datag, datab, 6, contrast, brightness);
  }else if(maptype == PPGPLOT_INVERTED_HEAT3) {
    float datal[] = {0.0, 0.2, 0.35, 0.5, 0.8, 1.0};
    float datar[] = {0.0, 0.0, 1.0, 1.0, 1.0, 1.0};
    float datag[] = {0.0, 0.0, 0.0, 0.0, 1.0, 1.0};
    float datab[] = {0.0, 1.0, 1.0, 0.0, 0.0, 1.0};
    ppgctab(datal, datar, datag, datab, 6, contrast, brightness);
  }else if(maptype == PPGPLOT_HEAT4) {
    float datal[] = {0.0, 0.33, 0.75, 1.0};
    float datar[] = {1.0, 0.0, 1.0, 1.0};
    float datag[] = {1.0, 0.0, 0.0, 0.0};
    float datab[] = {1.0, 1.0, 1.0, 0.0};
    ppgctab(datal, datar, datag, datab, 4, contrast, brightness);
  }else if(maptype == PPGPLOT_INVERTED_HEAT4) {
    float datal[] = {0.0, 0.33, 0.75, 1.0};
    float datar[] = {1.0, 1.0, 0.0, 1.0};
    float datag[] = {0.0, 0.0, 0.0, 1.0};
    float datab[] = {0.0, 1.0, 1.0, 1.0};
    ppgctab(datal, datar, datag, datab, 4, contrast, brightness);
  }else if(maptype == PPGPLOT_DIVREDBLUE) {
    float datal[] = {0.0, 0.5, 1.0};
    float datar[] = {1.0, 1.0, 0.0};
    float datag[] = {0.0, 1.0, 0.0};
    float datab[] = {0.0, 1.0, 1.0};
    ppgctab(datal, datar, datag, datab, 3, contrast, brightness);
  }else if(maptype == PPGPLOT_INVERTED_DIVREDBLUE) {
    float datal[] = {0.0, 0.5, 1.0};
    float datar[] = {0.0, 1.0, 1.0};
    float datag[] = {0.0, 1.0, 0.0};
    float datab[] = {1.0, 1.0, 0.0};
    ppgctab(datal, datar, datag, datab, 3, contrast, brightness);
  }else if(maptype == PPGPLOT_INFERNO) {
    float datal[] = {0.000000, 0.003922, 0.007843, 0.011765, 0.015686, 0.019608, 0.023529, 0.027451,
		     0.031373, 0.035294, 0.039216, 0.043137, 0.047059, 0.050980, 0.054902, 0.058824,
		     0.062745, 0.066667, 0.070588, 0.074510, 0.078431, 0.082353, 0.086275, 0.090196,
		     0.094118, 0.098039, 0.101961, 0.105882, 0.109804, 0.113725, 0.117647, 0.121569,
		     0.125490, 0.129412, 0.133333, 0.137255, 0.141176, 0.145098, 0.149020, 0.152941,
		     0.156863, 0.160784, 0.164706, 0.168627, 0.172549, 0.176471, 0.180392, 0.184314,
		     0.188235, 0.192157, 0.196078, 0.200000, 0.203922, 0.207843, 0.211765, 0.215686,
		     0.219608, 0.223529, 0.227451, 0.231373, 0.235294, 0.239216, 0.243137, 0.247059,
		     0.250980, 0.254902, 0.258824, 0.262745, 0.266667, 0.270588, 0.274510, 0.278431,
		     0.282353, 0.286275, 0.290196, 0.294118, 0.298039, 0.301961, 0.305882, 0.309804,
		     0.313725, 0.317647, 0.321569, 0.325490, 0.329412, 0.333333, 0.337255, 0.341176,
		     0.345098, 0.349020, 0.352941, 0.356863, 0.360784, 0.364706, 0.368627, 0.372549,
		     0.376471, 0.380392, 0.384314, 0.388235, 0.392157, 0.396078, 0.400000, 0.403922,
		     0.407843, 0.411765, 0.415686, 0.419608, 0.423529, 0.427451, 0.431373, 0.435294,
		     0.439216, 0.443137, 0.447059, 0.450980, 0.454902, 0.458824, 0.462745, 0.466667,
		     0.470588, 0.474510, 0.478431, 0.482353, 0.486275, 0.490196, 0.494118, 0.498039,
		     0.501961, 0.505882, 0.509804, 0.513725, 0.517647, 0.521569, 0.525490, 0.529412,
		     0.533333, 0.537255, 0.541176, 0.545098, 0.549020, 0.552941, 0.556863, 0.560784,
		     0.564706, 0.568627, 0.572549, 0.576471, 0.580392, 0.584314, 0.588235, 0.592157,
		     0.596078, 0.600000, 0.603922, 0.607843, 0.611765, 0.615686, 0.619608, 0.623529,
		     0.627451, 0.631373, 0.635294, 0.639216, 0.643137, 0.647059, 0.650980, 0.654902,
		     0.658824, 0.662745, 0.666667, 0.670588, 0.674510, 0.678431, 0.682353, 0.686275,
		     0.690196, 0.694118, 0.698039, 0.701961, 0.705882, 0.709804, 0.713725, 0.717647,
		     0.721569, 0.725490, 0.729412, 0.733333, 0.737255, 0.741176, 0.745098, 0.749020,
		     0.752941, 0.756863, 0.760784, 0.764706, 0.768627, 0.772549, 0.776471, 0.780392,
		     0.784314, 0.788235, 0.792157, 0.796078, 0.800000, 0.803922, 0.807843, 0.811765,
		     0.815686, 0.819608, 0.823529, 0.827451, 0.831373, 0.835294, 0.839216, 0.843137,
		     0.847059, 0.850980, 0.854902, 0.858824, 0.862745, 0.866667, 0.870588, 0.874510,
		     0.878431, 0.882353, 0.886275, 0.890196, 0.894118, 0.898039, 0.901961, 0.905882,
		     0.909804, 0.913725, 0.917647, 0.921569, 0.925490, 0.929412, 0.933333, 0.937255,
		     0.941176, 0.945098, 0.949020, 0.952941, 0.956863, 0.960784, 0.964706, 0.968627,
		     0.972549, 0.976471, 0.980392, 0.984314, 0.988235, 0.992157, 0.996078, 1.000000};
    float datar[] = {
      0.988362,0.982257,0.976511,0.971162,0.966249,0.961812,0.957896,0.954529,
      0.951740,0.949545,0.947937,0.946903,0.946403,0.946392,0.946809,0.947594,
      0.948683,0.950018,0.951546,0.953215,0.954997,0.956834,0.958720,0.960626,
      0.962517,0.964394,0.966243,0.968041,0.969783,0.971468,0.973088,0.974638,
      0.976108,0.977497,0.978806,0.980032,0.981173,0.982228,0.983196,0.984075,
      0.984865,0.985566,0.986175,0.986694,0.987124,0.987464,0.987714,0.987874,
      0.987945,0.987926,0.987819,0.987622,0.987337,0.986964,0.986502,0.985952,
      0.985315,0.984591,0.983779,0.982881,0.981895,0.980824,0.979666,0.978422,
      0.977092,0.975677,0.974176,0.972590,0.970919,0.969163,0.967322,0.965397,
      0.963387,0.961293,0.959114,0.956852,0.954506,0.952075,0.949562,0.946965,
      0.944285,0.941521,0.938675,0.935747,0.932737,0.929644,0.926470,0.923215,
      0.919879,0.916462,0.912966,0.909390,0.905735,0.902003,0.898192,0.894305,
      0.890341,0.886302,0.882188,0.878001,0.873741,0.869409,0.865006,0.860533,
      0.855992,0.851384,0.846709,0.841969,0.837165,0.832299,0.827372,0.822386,
      0.817341,0.812239,0.807082,0.801871,0.796607,0.791293,0.785929,0.780517,
      0.775059,0.769556,0.764010,0.758422,0.752794,0.747127,0.741423,0.735683,
      0.729909,0.724103,0.718264,0.712396,0.706500,0.700576,0.694627,0.688653,
      0.682656,0.676638,0.670599,0.664540,0.658463,0.652369,0.646260,0.640135,
      0.633998,0.627847,0.621685,0.615513,0.609330,0.603139,0.596940,0.590734,
      0.584521,0.578304,0.572081,0.565854,0.559624,0.553392,0.547157,0.540920,
      0.534683,0.528444,0.522206,0.515967,0.509730,0.503493,0.497257,0.491022,
      0.484789,0.478558,0.472328,0.466100,0.459875,0.453651,0.447428,0.441207,
      0.434987,0.428768,0.422549,0.416331,0.410113,0.403894,0.397674,0.391453,
      0.385228,0.379001,0.372768,0.366529,0.360284,0.354032,0.347771,0.341500,
      0.335217,0.328921,0.322610,0.316282,0.309935,0.303568,0.297178,0.290763,
      0.284321,0.277850,0.271347,0.264810,0.258234,0.251620,0.244967,0.238273,
      0.231538,0.224763,0.217949,0.211095,0.204209,0.197297,0.190367,0.183429,
      0.176493,0.169575,0.162689,0.155850,0.149073,0.142378,0.135778,0.129285,
      0.122908,0.116656,0.110536,0.104551,0.098702,0.092990,0.087411,0.081962,
      0.076637,0.071429,0.066331,0.061340,0.056449,0.051644,0.046915,0.042253,
      0.037668,0.033385,0.029432,0.025793,0.022447,0.019373,0.016561,0.013995,
      0.011663,0.009561,0.007676,0.006006,0.004547,0.003299,0.002267,0.001462};
    float datag[] = {
      0.998364,0.994109,0.989753,0.985282,0.980678,0.975924,0.971003,0.965896,
      0.960587,0.955063,0.949318,0.943348,0.937159,0.930761,0.924168,0.917399,
      0.910473,0.903409,0.896226,0.888942,0.881569,0.874129,0.866624,0.859069,
      0.851476,0.843848,0.836191,0.828515,0.820825,0.813122,0.805409,0.797692,
      0.789974,0.782258,0.774545,0.766837,0.759135,0.751442,0.743758,0.736087,
      0.728427,0.720782,0.713153,0.705540,0.697944,0.690366,0.682807,0.675267,
      0.667748,0.660250,0.652773,0.645320,0.637890,0.630485,0.623105,0.615750,
      0.608422,0.601122,0.593849,0.586606,0.579392,0.572209,0.565057,0.557937,
      0.550850,0.543798,0.536780,0.529798,0.522853,0.515946,0.509078,0.502249,
      0.495462,0.488716,0.482014,0.475356,0.468744,0.462178,0.455660,0.449191,
      0.442772,0.436405,0.430091,0.423831,0.417627,0.411479,0.405389,0.399359,
      0.393389,0.387481,0.381636,0.375856,0.370140,0.364492,0.358911,0.353399,
      0.347957,0.342586,0.337287,0.332060,0.326906,0.321827,0.316822,0.311892,
      0.307038,0.302260,0.297559,0.292933,0.288385,0.283913,0.279517,0.275197,
      0.270954,0.266786,0.262692,0.258674,0.254728,0.250856,0.247056,0.243327,
      0.239667,0.236077,0.232554,0.229097,0.225706,0.222378,0.219112,0.215906,
      0.212759,0.209670,0.206636,0.203656,0.200728,0.197851,0.195021,0.192239,
      0.189501,0.186807,0.184153,0.181539,0.178962,0.176421,0.173914,0.171438,
      0.168992,0.166575,0.164184,0.161817,0.159474,0.157151,0.154848,0.152563,
      0.150294,0.148039,0.145797,0.143567,0.141346,0.139134,0.136929,0.134729,
      0.132534,0.130341,0.128150,0.125960,0.123769,0.121575,0.119379,0.117179,
      0.114974,0.112764,0.110547,0.108322,0.106089,0.103848,0.101597,0.099338,
      0.097069,0.094790,0.092501,0.090203,0.087896,0.085580,0.083257,0.080927,
      0.078591,0.076253,0.073915,0.071579,0.069247,0.066925,0.064616,0.062325,
      0.060060,0.057827,0.055634,0.053490,0.051407,0.049396,0.047470,0.045644,
      0.043933,0.042353,0.040922,0.039647,0.038571,0.037705,0.037055,0.036621,
      0.036405,0.036405,0.036615,0.037030,0.037632,0.038400,0.039309,0.040329,
      0.041402,0.042489,0.043554,0.044559,0.045468,0.046242,0.046856,0.047293,
      0.047536,0.047574,0.047399,0.047008,0.046402,0.045583,0.044556,0.043328,
      0.041905,0.040294,0.038504,0.036590,0.034569,0.032474,0.030324,0.028139,
      0.025921,0.023702,0.021503,0.019331,0.017199,0.015133,0.013136,0.011225,
      0.009417,0.007713,0.006136,0.004692,0.003392,0.002249,0.001270,0.000466};
    float datab[] = {
      0.644924,0.631017,0.616760,0.602154,0.587206,0.571925,0.556275,0.540361,
      0.524203,0.507860,0.491426,0.474970,0.458592,0.442367,0.426373,0.410665,
      0.395289,0.380271,0.365627,0.351369,0.337475,0.323974,0.310820,0.298010,
      0.285546,0.273391,0.261534,0.249972,0.238686,0.227658,0.216877,0.206332,
      0.196018,0.185923,0.176037,0.166353,0.156863,0.147565,0.138453,0.129527,
      0.120785,0.112229,0.103863,0.095694,0.087731,0.079990,0.072489,0.065257,
      0.058329,0.051750,0.045581,0.039886,0.034916,0.030908,0.027814,0.025592,
      0.024202,0.023606,0.023770,0.024661,0.026250,0.028508,0.031409,0.034931,
      0.039050,0.043618,0.048392,0.053324,0.058367,0.063488,0.068659,0.073859,
      0.079073,0.084289,0.089499,0.094695,0.099874,0.105031,0.110164,0.115272,
      0.120354,0.125409,0.130438,0.135440,0.140417,0.145367,0.150292,0.155193,
      0.160070,0.164924,0.169755,0.174563,0.179350,0.184116,0.188860,0.193584,
      0.198286,0.202968,0.207628,0.212268,0.216886,0.221482,0.226055,0.230606,
      0.235133,0.239636,0.244113,0.248564,0.252988,0.257383,0.261750,0.266085,
      0.270390,0.274661,0.278898,0.283099,0.287264,0.291390,0.295477,0.299523,
      0.303526,0.307485,0.311399,0.315266,0.319085,0.322856,0.326576,0.330245,
      0.333861,0.337424,0.340931,0.344383,0.347777,0.351113,0.354388,0.357603,
      0.360757,0.363849,0.366879,0.369846,0.372748,0.375586,0.378359,0.381065,
      0.383704,0.386276,0.388781,0.391219,0.393589,0.395891,0.398125,0.400290,
      0.402385,0.404411,0.406369,0.408258,0.410078,0.411829,0.413511,0.415123,
      0.416667,0.418142,0.419549,0.420887,0.422156,0.423356,0.424488,0.425552,
      0.426548,0.427475,0.428334,0.429125,0.429846,0.430498,0.431080,0.431594,
      0.432039,0.432412,0.432714,0.432943,0.433098,0.433179,0.433183,0.433109,
      0.432955,0.432719,0.432400,0.431994,0.431497,0.430906,0.430217,0.429425,
      0.428524,0.427511,0.426377,0.425116,0.423721,0.422182,0.420491,0.418637,
      0.416608,0.414392,0.411976,0.409345,0.406485,0.403378,0.400007,0.396353,
      0.392400,0.388129,0.383522,0.378563,0.373238,0.367535,0.361447,0.354971,
      0.348111,0.340874,0.333277,0.325338,0.317085,0.308553,0.299776,0.290788,
      0.281624,0.272321,0.262912,0.253430,0.243904,0.234358,0.224813,0.215289,
      0.205799,0.196354,0.186962,0.177642,0.168414,0.159254,0.150164,0.141141,
      0.132232,0.123397,0.114621,0.105930,0.097327,0.088767,0.080282,0.071862,
      0.063460,0.055143,0.046836,0.038558,0.030909,0.024239,0.018570,0.013866};
    ppgctab(datal, datar, datag, datab, 256, contrast, brightness);
 }else {
    fflush(stdout);
    printerror(verbose.debug, "ERROR pgplot_set_maptype: Maptype undefined"); 
    return 0;
  }
  return 1;
}

//START REGION DEVELOP
//START REGION RELEASE

void clear_pgplot_frame(pgplot_frame_def_internal *frame)
{
  frame->svp  = 0;
  frame->svp_x1  = 0;
  frame->svp_x2  = 0;
  frame->svp_y1  = 0;
  frame->svp_y2  = 0;
  frame->swin = 0;
  frame->swin_showtwice = 0;
  frame->swin_x1  = 0;
  frame->swin_x2  = 0;
  frame->swin_y1  = 0;
  frame->swin_y2  = 0;
  frame->TR[0]  = 0;
  frame->TR[1]  = 0;
  frame->TR[2]  = 0;
  frame->TR[3]  = 0;
  frame->TR[4]  = 0;
  frame->TR[5]  = 0;
}

//START REGION DEVELOP
//START REGION RELEASE

/* Sets box options of pgplot plot to default values. */
void clear_pgplot_box(pgplot_box_definition *box)
{
  box->drawbox = 1;
  box->box_lw = 1;
  box->box_f = 1;
  box->box_labelsize = 1;
  box->box_xtick = 0; 
  box->box_ytick = 0;
  box->box_nxsub= 0;
  box->box_nysub= 0;
  sprintf(box->box_xopt, "bcnsti");
  sprintf(box->box_yopt, "bcnsti");
  box->drawtitle = 1;
  box->title_ch = 1;
  box->title_lw = 1;
  box->title_f = 1;
  box->title[0] = 0;
  box->drawlabels = 1;
  box->label_ch = 1;
  box->label_lw = 1;
  box->label_f = 1;
  box->dxlabel = 0;
  box->dylabel = 0;
  box->xlabel[0] = 0;
  box->ylabel[0] = 0;
  box->wedgelabel[0] = 0;
}

//START REGION DEVELOP
//START REGION RELEASE

/* Draw a box and labels.  */
void pgplot_drawbox(pgplot_box_definition *box)
{
  if(box->drawbox) {
    ppgslw(box->box_lw);
    ppgsch(box->box_labelsize);
    ppgscf(box->box_f);
    ppgbox(box->box_xopt,box->box_xtick,box->box_nxsub,box->box_yopt,box->box_ytick,box->box_nysub);
  }
  if(box->drawtitle) {
    ppgsch(box->title_ch);
    ppgslw(box->title_lw);
    ppgscf(box->title_f);
    
    ppgmtxt("T", (1.5*box->box_labelsize)/box->title_ch, 0.5, 0.5, box->title);
  }

  if(box->drawlabels) {
    ppgslw(box->label_lw);
    ppgsch(box->label_ch);
    ppgscf(box->label_f);
    /* Displacements is in units of the current character height.
       2.5: 1X for tick marks, 1X for actual numbers, 0.5 for extra margin */
    ppgmtxt("L", ((2.5+box->dylabel)*box->box_labelsize)/box->label_ch, 0.5, 0.5, box->ylabel);
    /* Displacements is in units of the current character height.
       1.5: 1X for actual numbers, 0.5 for extra margin 
        +1: displacement is for bottom of text, so have to displace with an extra character height
    */
    ppgmtxt("B", ((1.5+box->dxlabel)*box->box_labelsize)/box->label_ch+1, 0.5, 0.5, box->xlabel);
  }

  ppgsch(1);
  ppgslw(1);
  ppgscf(1);
}

//START REGION DEVELOP
//START REGION RELEASE

/*
  Plot 2D data (cmap) with nrx by nry points. viewport sets the size
  of the device and the placement of the graph within the device. The
  data should be uniformely sampled between xmin and xmax and ymin and
  ymax. You can zoom in by choosing xmingshow, xmaxshow
  differently. You can plot contours as well. In that case you can
  either provide nrcontours levels via contours, or set contours=NULL,
  in which case the contours are equally separated across the range of
  the data. When plotting contours, you might want to use the nogray
  option to only show the contours. The image transfer function can be
  set with itf (0=linear, 1=logarithmic, 2=square-root). You can use
  levelInversion to invert the colors (usefull in combination with
  itf). If forceMinZero is set than the black corresponds to 0, or
  else black is set to the minimum in the data. White is set to the
  maximum of the data if saturize is 1. Larger values of saturize mean
  that the peaks become clipped. Instead, if levelset is set, you can
  set the levels by levelmin and levelmax. The look of the labels etc
  is set with pgplotbox. If onlyData is set no labels/boxes are drawn
  (if onlyData is set to 2, the side panel boxes are drawn and the
  title if set). If sideright is set a collapsed panel is shown at the
  right. If set to 2, the numbers are not plotted. If set to 3 or 4,
  it is the same as 1 and 2, except that no tickmarkes are placed on
  axis along the map. If forceMinZeroRight is set the minimum of the
  scale in the side panel is set to zero. sidelw sets the linewidth
  used to draw the curve in the side panels. Showwedge indicates if
  there should be a color-range bar in the plot. If plotSubset is set,
  only the relevant bit of the map will be sent to pgplot, thereby
  making ps files potentially a lot smaller. If debug is set info
  about the color scale is printed to the stdout. The xmin,xmax and
  ymin,ymax corresponds to the centres of the relevant bins. If
  showTwice is set the map is plotted twice above each other.

  The color maps are:
    PPGPLOT_GRAYSCALE
    PPGPLOT_INVERTED_GRAYSCALE
    etc.

  Return 1: Succesful
         0: Plot device couldn't be opened or no memory could be allocated
*/

int pgplotMap(pgplot_options_definition *pgplot, float *cmap, int nrx, int nry, float xmin, float xmax, float xminshow, float xmaxshow, float ymin, float ymax, float yminshow, float ymaxshow, int maptype, int itf, int nogray, int nrcontours, float *contours, int contourlw, int forceMinZero, float saturize, int levelset, float levelmin, float levelmax, int levelInversion, int onlyData, int sideright, int forceMinZeroRight, int sidetop, int forceMinZeroTop, int sidelw, int showwedge, int plotSubset, int showTwice, verbose_definition verbose)
{
  float datamin, datamax, xmin2, xmax2;
  float *collapse, x, y, remember2_x1, remember2_x2, remember2_y1, remember2_y2;
  int i, j, deviceID, firstc, lastc, firstpoint;
  float junk_f, lasty;
  pgplot_frame_def_internal pgplot_frame_internal;
  pgplot_options_definition *pgplot_backup;

  pgplot_backup = (pgplot_options_definition *)malloc(sizeof(pgplot_options_definition));
  if(pgplot_backup == NULL) {
    printerror(verbose.debug, "ERROR pgplotMap: Memory allocation error");
    return 0;
  }
  memcpy(pgplot_backup, pgplot, sizeof(pgplot_options_definition));

  clear_pgplot_frame(&pgplot_frame_internal);

  /* Remember mapping for pgplotMapCoordinate and pgplot_makeframe */
  internal_pgplot_xmin = xmin;
  internal_pgplot_xmax = xmax;
  internal_pgplot_ymin = ymin;
  internal_pgplot_ymax = ymax;
  internal_pgplot_nrx = nrx;
  internal_pgplot_nry = nry;

  // Not necessary, but gets rid of warning when optimising code
  lasty = 0;

  if(verbose.debug) {
    printf("pgplotMap -- dimensions: %dX%d points  xrange=%e %e yrange=%e %e\n", nrx, nry, xmin, xmax, ymin, ymax);
    printf("pgplotMap -- Horizontal dimensions %f - %f (shown %f - %f)\n", xmin, xmax, xminshow, xmaxshow);
    printf("pgplotMap -- Vertical dimensions %f - %f (shown %f - %f)\n", ymin, ymax, yminshow, ymaxshow);
  }
  /* Open device, unless want to re-use already opened device. */
  if(pgplot_opendevice(&(pgplot->viewport), &deviceID, verbose) == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR pgplotMap: Cannot open plot device");
    memcpy(pgplot, pgplot_backup, sizeof(pgplot_options_definition));
    free(pgplot_backup);
    return 0;
  }
  if(onlyData == 0) {
    pgplot_setWindowsize(pgplot->viewport.windowwidth, pgplot->viewport.windowheight, pgplot->viewport.aspectratio);
  }

  /* Reset some pgplot settings */
  ppgslw(1);
  ppgsci(1);
  ppgscf(1);
  ppgsch(1);

  /* Calculate the viewport coordinates for the main panel if required */
  if(onlyData == 0 || onlyData == 2) {
    /*
    if(sideright) 
      pgplot_frame_internal.svp_x2 = (0.75-pgplot->viewport.dxplot-0.15)*pgplot->viewport.xsize+pgplot->viewport.dxplot+0.15;
    else
      pgplot_frame_internal.svp_x2 = (0.90-pgplot->viewport.dxplot-0.15)*pgplot->viewport.xsize+pgplot->viewport.dxplot+0.15;
    if(sidetop) 
      pgplot_frame_internal.svp_y2 = (0.75-pgplot->viewport.dyplot-0.15)*pgplot->viewport.ysize+pgplot->viewport.dyplot+0.15;
    else
      pgplot_frame_internal.svp_y2 = (0.90-pgplot->viewport.dyplot-0.15)*pgplot->viewport.ysize+pgplot->viewport.dyplot+0.15;
    pgplot_frame_internal.svp_x1 = pgplot->viewport.dxplot+0.15 + pgplot->viewport.dxplot;
    pgplot_frame_internal.svp_y1 = pgplot->viewport.dyplot+0.15 + pgplot->viewport.dxplot;
    */
    pgplot_frame_internal.svp_x1 = pgplot->viewport.dxplot+0.15;
    pgplot_frame_internal.svp_y1 = pgplot->viewport.dyplot+0.15;
    if(sideright) 
      pgplot_frame_internal.svp_x2 = (0.75-0.15)*pgplot->viewport.xsize+pgplot_frame_internal.svp_x1;
    else
      pgplot_frame_internal.svp_x2 = (0.90-0.15)*pgplot->viewport.xsize+pgplot_frame_internal.svp_x1;
    if(sidetop) 
      pgplot_frame_internal.svp_y2 = (0.75-0.15)*pgplot->viewport.ysize+pgplot_frame_internal.svp_y1;
    else
      pgplot_frame_internal.svp_y2 = (0.90-0.15)*pgplot->viewport.ysize+pgplot_frame_internal.svp_y1;

    //    printf("YYYYY %f %f %f %f\n", pgplot_frame_internal.svp_x1, pgplot_frame_internal.svp_x2, pgplot_frame_internal.svp_y1, pgplot_frame_internal.svp_y2);
    pgplot_frame_internal.svp = 1;
  }else {
    pgplot_frame_internal.svp = 0;
  }
  /* Set the coordinate ranges of the axes. */
  pgplot_frame_internal.swin = 1;
  pgplot_frame_internal.swin_showtwice = showTwice;
  pgplot_frame_internal.swin_x1 = xminshow;
  pgplot_frame_internal.swin_x2 = xmaxshow;
  pgplot_frame_internal.swin_y1 = yminshow;
  pgplot_frame_internal.swin_y2 = ymaxshow;
  if(showTwice)
    pgplot_frame_internal.swin_y2 += ymaxshow-yminshow;
  /* Create the main panel */
  pgplot_makeframe(&pgplot_frame_internal);


  /* Calculate range in data values within the confinement of the axis. */
  firstpoint = 1;
  datamin = datamax = 0;
  for(i = 0; i < nrx; i++) {
    for(j = 0; j <nry; j++) {
      /* Calculate the axis coordinates of the point */
      pgplotMapCoordinateInverse(&x, &y, i, j);
      //      printf("XXXXXX x=%f\n", x);
      if(x >= xminshow && x <= xmaxshow && y >= yminshow && y <= ymaxshow) {
	if(firstpoint) {
	  datamin = datamax = cmap[j*nrx+i];
	  firstpoint = 0;
	}
	if(cmap[j*nrx+i] > datamax)
	  datamax = cmap[j*nrx+i];
	if(cmap[j*nrx+i] < datamin) {
	  datamin = cmap[j*nrx+i];
	}
      }
    }
  }
  if(verbose.debug) printf("pgplotMap -- Data range: %e to %e\n", datamin, datamax);
  if(forceMinZero)
    datamin = 0;
  datamax = datamin + (datamax-datamin)/saturize;
  if(levelset) {
    datamax = levelmax;
    datamin = levelmin;
  }
  if(levelInversion) {
    junk_f = datamax;
    datamax = datamin;
    datamin = junk_f;
  }
  //  fprintf(stderr, "levelset=%d, levelInversion=%d\n", levelset, levelInversion);
  if(verbose.debug) printf("pgplotMap -- Color range: %e to %e\n", datamin, datamax);


  /* Find out what range of color indices is supported by device */
  firstc = 17;
  lastc = 1000;
  ppgscir(firstc, lastc);
  cpgqcir(&firstc, &lastc);
  if(verbose.debug) printf("pgplotMap -- Color index range:  %d to %d\n", firstc, lastc);

  /* Set the image transfer function */
  ppgsitf(itf);

  if(pgplot_set_maptype(maptype, verbose) == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR pgplotMap: Cannot set colour map"); 
    memcpy(pgplot, pgplot_backup, sizeof(pgplot_options_definition));
    free(pgplot_backup);
    return 0;
  }


  /* Find the subset of the map which is actually plotted. This can
     drastically reduce the file size of the generated postscript
     file. */
  int subset_nrx, subset_nry, subset_extrabefore, subset_extraafter;
  float *cmap_subset;


  if(plotSubset) {
    int subset_x0, subset_y0;
    //    float subset_xmin, subset_xmax, subset_ymin, subset_ymax
    
    /* Work out some numbers that characterise the new dimensions of the smaller map */
    subset_extrabefore = 1;
    subset_extraafter = 1;
    subset_x0 = (nrx-1)*(xminshow - xmin)/(xmax-xmin)-1;
    if(subset_x0 < 0)
      subset_x0 = 0;
    //    printf("XXXXXXX nrx=%d xmaxshow=%f xmax=%f xmin=%f subset_x0=%d\n", nrx, xmaxshow, xmax, xmin, subset_x0);
    if(nrx != 1) {
      subset_nrx = (nrx-1)*(xmaxshow - xmin)/(xmax-xmin)+1 - subset_x0;
    }else {
      subset_nrx = 1 - subset_x0;
    }
    if(subset_nrx > nrx)
      subset_nrx = nrx;
    
    subset_y0 = (nry-1)*(yminshow - ymin)/(ymax-ymin)-1;
    if(subset_y0 < 0) {
      subset_y0 = 0;
      subset_extrabefore = 0;
    }
    if(nry != 1)
      subset_nry = (nry-1)*(ymaxshow - ymin)/(ymax-ymin)+1 - subset_y0;
    else
      subset_nry = +1 - subset_y0; 
    if(subset_nry >= nry) {
      subset_nry = nry;
      subset_extraafter = 0;
    }
      
    /* Check if we really want to generate a subset */
    if(subset_nrx == nrx && subset_nry == nry && subset_y0 == 0 && subset_x0 == 0)
      plotSubset = 0;
    if(verbose.debug && plotSubset) printf("pgplotMap -- changed to: %dX%d points  startx=%d starty=%d\n", subset_nrx, subset_nry, subset_x0, subset_y0);
    
    /* Only generate subset if it is required */
    if(plotSubset) {
      /* Try to allocate memory for the new subset map */
      cmap_subset = (float *)malloc(subset_nrx*subset_nry*sizeof(float));
      if(cmap_subset == NULL) {
	fflush(stdout);
	printwarning(verbose.debug, "WARNING pgplotMap: Cannot allocate memory to plot subset of the map");
	if(verbose.debug) {
	  printwarning(verbose.debug, "Tried to allocate %ld x %ld floating points", subset_nrx, subset_nry);
	}
	plotSubset = 0;
      }
    }
    
    /* Only generate subset if it is required and if memory could be allocated */
    if(plotSubset) {
      /* Generate subset */
      for(i = 0; i < subset_nrx; i++) {
	for(j = 0; j < subset_nry; j++) {
	  cmap_subset[j*subset_nrx+i] = cmap[(j+subset_y0)*nrx+i+subset_x0];
	}
      }
      /*      
	      subset_xmin = subset_x0*(xmax-xmin)/(float)(nrx-1) + xmin;
	      subset_xmax = (subset_x0+subset_nrx)*(xmax-xmin)/(float)(nrx-1) + xmin;
	      subset_ymin = subset_y0*(ymax-ymin)/(float)(nry-1) + ymin;
	      subset_ymax = (subset_y0+subset_nry)*(ymax-ymin)/(float)(nry-1) + ymin;
	      printf("XXXX range: %f to %f and %f to %f\n", subset_xmin, subset_xmax, subset_ymin, subset_ymax); 
      */
      pgplot_frame_internal.TR[0] += pgplot_frame_internal.TR[1]*subset_x0;
      pgplot_frame_internal.TR[3] += pgplot_frame_internal.TR[5]*subset_y0;
      /*      TR[0] = subset_xmin - (subset_xmax-subset_xmin)/(float)(subset_nrx-1); 
	      TR[1] = (subset_xmax-subset_xmin)/(float)(subset_nrx-1);    TR[2] = 0;
	      TR[3] = subset_ymin - (subset_ymax-subset_ymin)/(float)(subset_nry-1); TR[4] = 0;    TR[5] = (subset_ymax-subset_ymin)/(float)(subset_nry-1); */
      /*      printf("XXXX TR: %f %f %f %f %f %f\n", TR[0], TR[1], TR[2], TR[3], TR[4], TR[5]); */
    }
  }

  /* If we didn't make a subset, make the subset identical to the original */
  if(plotSubset == 0) {
    subset_nrx = nrx;
    subset_nry = nry;
    subset_extrabefore = 0;
    subset_extraafter = 0;
    cmap_subset = cmap;
  }

  /* Draw the data map */
  if(!nogray) {
    if(showTwice) {
      pgplot_frame_internal.TR[3] += pgplot_frame_internal.TR[5]*(subset_nry-(subset_extrabefore + subset_extraafter));
      ppgimag(cmap_subset, subset_nrx, subset_nry, 1, subset_nrx, 1, subset_nry, datamax, datamin, pgplot_frame_internal.TR);
      pgplot_frame_internal.TR[3] -= pgplot_frame_internal.TR[5]*(subset_nry-(subset_extrabefore + subset_extraafter));
    }
    //    printf("XXXXXX PLOTTING %p %p\n", cmap_subset, cmap);
    //    printf("XXXXXX nrx/y = %d %d\n", subset_nrx, subset_nry);
    //    printf("XXXXXX datamax/min = %e %e\n", datamin, datamax);
    //    printf("XXXXXX TR=%e %e %e %e %e %e\n", pgplot_frame_internal.TR[0], pgplot_frame_internal.TR[1], pgplot_frame_internal.TR[2], pgplot_frame_internal.TR[3], pgplot_frame_internal.TR[4], pgplot_frame_internal.TR[5]);
    //    printf("XXXXXX Values are %e %e %e %e.....\n", cmap_subset[0], cmap_subset[1], cmap_subset[2], cmap_subset[3]);
    ppgimag(cmap_subset, subset_nrx, subset_nry, 1, subset_nrx, 1, subset_nry, datamax, datamin, pgplot_frame_internal.TR);
  }
  /* Create a wedge showing the color mapping */
  if(showwedge && sideright == 0) {
    ppgsch(pgplot->box.box_labelsize*0.5);
    ppgslw(ceil(0.5*pgplot->box.box_lw));
    ppgwedg("RI", 0, 6, datamax, datamin, pgplot->box.wedgelabel);
    ppgslw(1);
  }
  /* Show contours */
  if(nrcontours > 0) {
    int didallocate_contours = 0;
    if(contours == NULL) {   // User didn't provide the contour levels, so set them ourselves
      contours = malloc(nrcontours*sizeof(float));
      if(contours == NULL) {
	printerror(verbose.debug, "ERROR pgplotMap: Cannot allocate memory");
	return 0;
      }
      didallocate_contours = 1;
      int contlevel;
      for(contlevel = 0; contlevel < nrcontours; contlevel++) {
	contours[contlevel] = datamin + (contlevel+1)*(datamax-datamin)/(float)(nrcontours+1);
      }
    }
    ppgslw(contourlw);
    if(showTwice) {
      pgplot_frame_internal.TR[3] += pgplot_frame_internal.TR[5]*(subset_nry-(subset_extrabefore + subset_extraafter));
      ppgcont(cmap_subset, subset_nrx, subset_nry, 1, subset_nrx, 1, subset_nry, contours, nrcontours, pgplot_frame_internal.TR);
      pgplot_frame_internal.TR[3] -= pgplot_frame_internal.TR[5]*(subset_nry-(subset_extrabefore + subset_extraafter));
    }
    ppgcont(cmap_subset, subset_nrx, subset_nry, 1, subset_nrx, 1, subset_nry, contours, nrcontours, pgplot_frame_internal.TR);
    if(didallocate_contours) {
      free(contours);
      contours = NULL;
    }
  }
  
  /* Clear memory of the subset if required */
  if(plotSubset)
    free(cmap_subset);


  /* Generate the labels etc. */
  if(onlyData != 0)
    pgplot->box.drawbox = 0;
  if(sidetop && onlyData != 1)
    pgplot->box.drawtitle = 0;
  if(onlyData != 0)
    pgplot->box.drawlabels = 0;

  pgplot_drawbox(&(pgplot->box));

  if(sideright) {
    collapse = (float *)calloc(nry, sizeof(float));
    if(collapse == NULL) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR pgplotMap: Cannot allocate memory");
      memcpy(pgplot, pgplot_backup, sizeof(pgplot_options_definition));
      free(pgplot_backup);
      return 0;
    }
    //    printf("XXXXX nry=%d\n", nry);
    for(i = 0; i < nrx; i++) {
      x = i*(xmax-xmin)/(float)(nrx-1) + xmin;
      if(x >= xminshow && x <= xmaxshow) {
	for(j = 0; j <nry; j++) {
	  collapse[j] += cmap[j*nrx+i];
	}
      }
    }
    //    printf("XXXXX yshow=%f %f\n", yminshow, ymaxshow);
    xmin2 = xmax2 = collapse[0];
    i = 0;
    for(j = 0; j <nry; j++) {
      y = j*(ymax-ymin)/(float)(nry-1) + ymin;
      if(y >= yminshow && y <= ymaxshow) {
	if(i == 0)
	  xmin2 = xmax2 = collapse[j];
	if(collapse[j] > xmax2)
	  xmax2 = collapse[j];
	if(collapse[j] < xmin2)
	  xmin2 = collapse[j];
	i++;
      }
    }
    //    printf("XXXXX sideright %f %f\n", xmin2, xmax2);
    if(xmin2 == xmax2) {
      float offset;
      offset = 0.05*xmax2;
      xmin2 -= offset;
      xmax2 += offset;
      if(xmin2 == xmax2) {
	xmin2 -= 0.5;
	xmax2 += 0.5;
      }
    }
    //    printf("XXXXX sideright %f %f (onlyData=%d sideright=%d)\n", xmin2, xmax2, onlyData, sideright);
    if(forceMinZeroRight)
      xmin2 = 0;
    if(sidetop) 
      junk_f = (0.75-0.15)*pgplot->viewport.ysize+pgplot->viewport.dyplot+0.15;
    else
      junk_f = (0.90-0.15)*pgplot->viewport.ysize+pgplot->viewport.dyplot+0.15;
    if(onlyData != 0) {
      ppgqvp(0, &pgplot_frame_internal.svp_x1, &pgplot_frame_internal.svp_x2, &pgplot_frame_internal.svp_y1, &pgplot_frame_internal.svp_y2);
      ppgqwin(&remember2_x1, &remember2_x2, &remember2_y1, &remember2_y2);
    }
    ppgsvp(pgplot_frame_internal.svp_x2, pgplot_frame_internal.svp_x2+0.15, pgplot_frame_internal.svp_y1, junk_f);
    if(showTwice == 0)
      ppgswin(xmin2, xmax2, yminshow, ymaxshow);
    else
      ppgswin(xmin2, xmax2, yminshow, ymaxshow+ymaxshow-yminshow);
    if(showwedge) {
      ppgsch(pgplot->box.box_labelsize*0.5);
      //      printf("XXXX: set lw to %d = half %d\n", ceil(0.5*pgplot->box.box_lw), pgplot->box.box_lw);
      ppgslw(ceil(0.5*pgplot->box.box_lw));
      ppgwedg("RI", 0, 6, datamax, datamin, pgplot->box.wedgelabel);
      ppgslw(1);
    }
    if(onlyData == 0 || onlyData == 2) {
      ppgslw(pgplot->box.box_lw);
      ppgsch(pgplot->box.box_labelsize*0.5);
      if(sideright == 4) {
	ppgbox("bcsti",0.0,0,"bc",0.0,0);
      }else if(sideright == 3) {
	ppgbox("bcnsti",0.0,0,"bc",0.0,0);
      }else if(sideright == 2) {
	ppgbox("bcsti",0.0,0,"bctsi",0.0,0);
      }else {
	ppgbox("bcnsti",0.0,0,"bctsi",0.0,0);
      }
      ppgsch(1);
      ppgslw(1);
    }
    ppgslw(sidelw);
    for(j = 0; j < nry; j++) {
      x = collapse[j];
      y = j*(ymax-ymin)/(float)(nry-1) + ymin;
      if(y >= yminshow && y <= ymaxshow) {
	if(j == 0)
	  ppgmove(x, y);
	else {
	  ppgdraw(x, y);
	  lasty = y;
	}
      }else {
	  ppgmove(x, y);
      }
      //      printf("XXXXX  DRAW %f %f (%f %f)\n", x, y, ymin, ymax);
    }
    if(showTwice) {
      for(j = 0; j < nry; j++) {
	x = collapse[j];
	if(plotSubset == 0)
	  y = j*(ymax-ymin)/(float)(nry-1) + ymin + nry*(ymax-ymin)/(float)(nry-1);
	else
	  y = j*(ymax-ymin)/(float)(nry-1) + ymin + pgplot_frame_internal.TR[5]*(subset_nry-(subset_extrabefore + subset_extraafter));
	if(y <= ymaxshow+ymaxshow-yminshow) {
	  if(y > lasty)
	    ppgdraw(x, y);
	}
	ppgsci(1);
      }
    }
    ppgslw(1);
    free(collapse);
    if(onlyData != 0) {
      ppgsvp(pgplot_frame_internal.svp_x1, pgplot_frame_internal.svp_x2, pgplot_frame_internal.svp_y1, pgplot_frame_internal.svp_y2);
      ppgswin(remember2_x1, remember2_x2, remember2_y1, remember2_y2);
    }
  }

  if(sidetop) {
    collapse = (float *)calloc(nrx, sizeof(float));
    if(collapse == NULL) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR pgplotMap: Cannot allocate memory");
      memcpy(pgplot, pgplot_backup, sizeof(pgplot_options_definition));
      free(pgplot_backup);
      return 0;
    }
    for(i = 0; i < nry; i++) {
      if(nry == 1)
	y = ymin;
      else
	y = i*(ymax-ymin)/(float)(nry-1) + ymin;
      if(y >= yminshow && y <= ymaxshow) {
	for(j = 0; j <nrx; j++) {
	  collapse[j] += cmap[i*nrx+j];
	}
      }
    }
    i = 0;
    xmin2 = xmax2 = 0;
    for(j = 0; j <nrx; j++) {
      x = (j*(xmax-xmin)/(float)(nrx - 1)) + xmin;
      if(x >= xminshow && x <= xmaxshow) {
	if(i == 0) {
	  xmin2 = xmax2 = collapse[j];
	  i = 1;
	}
	if(collapse[j] > xmax2)
	  xmax2 = collapse[j];
	if(collapse[j] < xmin2)
	  xmin2 = collapse[j];
      }
    }
    //    printf("XXXXX sidetop %f %f\n", xmin2, xmax2);
    if(xmin2 == xmax2) {
      float offset;
      offset = 0.05*xmax2;
      xmin2 -= offset;
      xmax2 += offset;
      if(xmin2 == xmax2) {
	xmin2 -= 0.5;
	xmax2 += 0.5;
      }
    }
    //    printf("XXXXX sidetop %f %f (onlyData=%d sidetop=%d)\n", xmin2, xmax2, onlyData, sideright);
    if(forceMinZeroRight)
      xmin2 = 0;
    if(sideright)
      junk_f = (0.75-0.15)*pgplot->viewport.xsize+pgplot->viewport.dxplot+0.15;
    else
      junk_f = (0.90-0.15)*pgplot->viewport.xsize+pgplot->viewport.dxplot+0.15;
    if(onlyData != 0) {
      ppgqvp(0, &pgplot_frame_internal.svp_x1, &pgplot_frame_internal.svp_x2, &pgplot_frame_internal.svp_y1, &pgplot_frame_internal.svp_y2);
      ppgqwin(&remember2_x1, &remember2_x2, &remember2_y1, &remember2_y2);
      junk_f = pgplot_frame_internal.svp_x2;
    }
    ppgsvp(pgplot_frame_internal.svp_x1, junk_f, pgplot_frame_internal.svp_y2, pgplot_frame_internal.svp_y2+0.15);
    ppgswin(xminshow, xmaxshow, xmin2, xmax2);
    if(onlyData == 0 || onlyData == 2) {
      ppgsch(pgplot->box.box_labelsize*0.5);
      ppgslw(pgplot->box.box_lw);
      if(sidetop == 2) {
	ppgbox("bcsti",0.0,0,"bctsi",0.0,0);
      }else {
	ppgbox("bcsti",0.0,0,"bcntsi",0.0,0);
      }
      ppgsch(1);
      ppgslw(1);
    }
    ppgslw(sidelw);
    for(j = 0; j < nrx; j++) {
      y = collapse[j];
      x = j*(xmax-xmin)/(float)(nrx-1) + xmin;
      if(j == 0)
	ppgmove(x, y);
      else
	ppgdraw(x, y);
    }
    ppgslw(1);
    if(onlyData == 0 || onlyData == 2) {
      /*
      printf("XXXX second %d %d\n", onlyData, sidetop);
      printf("XXXX yminshow=%f ymaxshow=%f labelsize=%f\n", yminshow, ymaxshow, labelsize);
      printf("XXXX pos=%f '%s'\n", ymaxshow+0.06*(labelsize/1.7)*(ymaxshow-yminshow), title);
      */
      ppgsch(pgplot->box.title_ch);
      ppgslw(pgplot->box.title_lw);
      ppgscf(pgplot->box.title_f);
      ppgptxt(0.5*(xmaxshow-xminshow)+xminshow, xmax2+0.3*(pgplot->box.label_ch/1.7)*(xmax2-xmin2), 0, 0.5, pgplot->box.title);
      ppgscf(1);
      ppgslw(1);
      ppgsch(pgplot->box.label_ch);
    }
    free(collapse);
    if(onlyData != 0) {
      ppgsvp(pgplot_frame_internal.svp_x1, pgplot_frame_internal.svp_x2, pgplot_frame_internal.svp_y1, pgplot_frame_internal.svp_y2);
      ppgswin(remember2_x1, remember2_x2, remember2_y1, remember2_y2);
    }
  }

  /* Reset the viewport to the main panel, which is useful for interactive programs */
  if(onlyData == 0 || onlyData == 2) {
    ppgsvp(pgplot_frame_internal.svp_x1, pgplot_frame_internal.svp_x2, pgplot_frame_internal.svp_y1, pgplot_frame_internal.svp_y2);
    if(showTwice == 0)
      ppgswin(xminshow, xmaxshow, yminshow, ymaxshow);
    else
      ppgswin(xminshow, xmaxshow, yminshow, ymaxshow+ymaxshow-yminshow);
  }

  memcpy(pgplot, pgplot_backup, sizeof(pgplot_options_definition));
  free(pgplot_backup);
  if(!pgplot->viewport.dontclose)
    ppgclos();
  return 1;
}

//START REGION DEVELOP
//START REGION RELEASE


/* Ask user to define pulse longitude ranges in a graph with nrBins
   bins. If onlyOne is specified the user can keep on selecting one
   component until he/she is satisfied. If powerTwo is selected the
   region is forced to be a power of two. If evenNumber is set, the
   region is forced to be an even number. In viewport, if dontclose is
   set the plotdevice is not closed, so you can do additional
   drawing. If dontopen is set it assumes a plotting device is opened
   and selected. It returns the number of selected regions, which can
   be zero if nothing is selected. -1 is returned if device doesn't
   support a cursor. */
int selectRegions(float *profileI, int nrBins, pgplot_options_definition *pgplot, int onlyOne, int powerTwo, int evenNumber, pulselongitude_regions_definition *regions, verbose_definition verbose)
{
  int i, j, k, bin1, bin2, replot;
  float x, y;
  char c;
  pgplot_options_definition *pgplot2;

  pgplot2 = (pgplot_options_definition *)malloc(sizeof(pgplot_options_definition));
  if(pgplot2 == NULL) {
    printerror(verbose.debug, "ERROR selectRegions: Memory allocation error");
    return 0;
  }
  x = y = 0;
  c = 0;

  // Just to get rid of compiler warning when optimising code
  bin1 = bin2 = 0;

  if(strcmp(pgplot->viewport.plotDevice, "?") == 0)
    printf("Select plotting device to show profile in order to select the onpulse region\n  ");
  memcpy(pgplot2, pgplot, sizeof(pgplot_options_definition));
  pgplot2->viewport.dontclose = 1; /* Never close device at this point */
  /*
  pgplot_box_def pgplotbox;
  clear_pgplot_box(&pgplotbox);
  strcpy(pgplot->box.xlabel, "Bin");
  strcpy(pgplot->box.ylabel, "Intensity");
  strcpy(pgplot->box.title, title);
  */

  pgplotGraph1(pgplot2, profileI, NULL, NULL, nrBins, 0, nrBins-1, 0, 0, nrBins-1, 0, 0, 0, 1, 0, 1, 0, 1, 1, regions, -1, verbose);
  printf("Select onpulse region (click left and right of region, you can select multiple regions and pressing 'r' inside the PGPLOT window removes selections one by one).\nWhen finished press S inside the PGPLOT window.\n");
  i = 0;
  replot = 0;
  do {
    if(i == 0)
      cpgband(6, 0, 0.0, 0.0, &x, &y, &c);
    else
      cpgband(4, 0, bin1-0.5, 0.0, &x, &y, &c); 
    if(c == 0) {
      fflush(stdout);
      printwarning(verbose.debug, "WARNING selectRegions: PGPLOT device does not support cursor: no pulse phase selection can be made.");
      if(pgplot->viewport.dontclose == 0) {
	ppgclos();
	fflush(stdout);
      }
      free(pgplot2);
      return -1;          // No cursor supported in device
    }else if(c == 65) {    /* left mouse button */
      if(i == 0)
	bin1 = (int)(x+0.5);
      else
	bin2 = (int)(x+0.5);
	//	bin2 = (int)ceil(x);
      i++;
    }else if(c == 115 || c == 83 ) {   /* s or S*/
      i = 3;
      printf("S is pressed\n");
    }else if(c == 'r' || c == 'R') {
      printf("R is pressed\n");
      i = 0;
      regions->nrRegions--;
      if(regions->nrRegions < 0)
	regions->nrRegions = 0;
      replot = 1;
    }
    if(i == 2) {             /* Both edges defined */
      if(bin2 < bin1) {
	i = bin1;
	bin1 = bin2;
	bin2 = i;
      }
      if(bin1 < 0)
	bin1 = 0;
      if(bin2 >= nrBins)
	bin2 = nrBins-1;
      if(verbose.verbose) printf("Selected region: %d - %d\n", bin1, bin2);

      if(powerTwo == 1) {    /* Make selected region power of 2 */
  	y = log(bin2-bin1+1.0)/log(2.0);
	j = y;
	k = pow(2.0, j+1.0);
	bin1 -= 0.5*(k - (bin2-bin1+1.0));
	if(bin1 < 0)
	  bin1 = 0;
	bin2 = bin1+k-1;
	if(verbose.verbose) printf("Selected region: %d - %d\n", bin1, bin2);
      }else if(evenNumber == 1) {
	j = bin2-bin1+1;
	j /= 2;
	j = 2*j-(bin2-bin1+1);
	if(j != 0) {   /* make even number */
	  bin2 += 1;
	}
	if(bin2 >= nrBins) {
	  bin2 = nrBins-1;
	  bin1 -= 1;
	}
	if(bin1 < 0)
	  bin1 = 0;
	if(verbose.verbose) printf("Selected region: %d - %d\n", bin1, bin2);
      }

      if(onlyOne)
	regions->nrRegions = 0;
      regions->left_bin[regions->nrRegions] = bin1;
      regions->right_bin[regions->nrRegions] = bin2;
      regions->bins_defined[regions->nrRegions] = 1;
      regions->nrRegions += 1;
      if(regions->nrRegions == MAX_pulselongitude_regions) {
	fflush(stdout);
	printerror(verbose.debug, "ERROR selectRegions: Too many regions selected.");
	ppgend();
	free(pgplot2);
	return regions->nrRegions;
      }
      replot = 1;
      i = 0;
    }
    if(replot) {
      memcpy(pgplot2, pgplot, sizeof(pgplot_options_definition));
      pgplot2->viewport.dontclose = 1; /* Never close device at this point */
      pgplot2->viewport.dontopen  = 1; /* Never open device at this point */
      pgplotGraph1(pgplot2, profileI, NULL, NULL, nrBins, 0, nrBins-1, 0, 0, nrBins-1, 0, 0, 0, 1, 0, 1, 0, 1, 1, regions, -1, verbose);
      replot = 0;
    }
  }while(i < 3);
  if(pgplot->viewport.dontclose == 0)
    ppgclos();
  free(pgplot2);
  return regions->nrRegions;
}

//START REGION DEVELOP
//START REGION RELEASE

/* Check if bin is in the region specified by whichregion (counting
   from 1). Zero is returned when bin is outside the specified
   whichregion. If whichregion is 0, then it will return the first
   region (counting from 1) in which bin falls (or 0 if it is in none
   of the regions). If regions set to NULL it will return 0.*/
int checkRegions(int bin, pulselongitude_regions_definition *regions, int whichregion, verbose_definition verbose)
{
  int i;
  if(regions == NULL)
    return 0;
  for(i = 0; i < regions->nrRegions; i++) {
    if(regions->bins_defined[i] == 0) {
      fflush(stdout);
      printwarning(verbose.debug, "WARNING: checkRegions called without region defined in bins, call is ignored");
      return 0;
    }
    if(bin >= regions->left_bin[i] && bin <= regions->right_bin[i]) {
      if(i+1 == whichregion || whichregion == 0) 
	 return i+1;
    }
  }
  return 0;
}

//START REGION DEVELOP

/* Prints the regions defined in regions. Se the nr of spaces before
   the output by indent. */
void printRegions(pulselongitude_regions_definition *regions, verbose_definition verbose)
{
  int n, i;
  if(regions == NULL) {
    for(i = 0; i < verbose.indent; i++)
      printf(" ");
    printf("Regions are not defined\n");
    return;
  }
  for(i = 0; i < verbose.indent; i++)
    printf(" ");
  printf("There are %d regions\n", regions->nrRegions);
  for(n = 0; n < regions->nrRegions; n++) {
    if(regions->bins_defined[n]) {
      for(i = 0; i < verbose.indent; i++)
	printf(" ");
      printf("  Region %d: left bin=%d, right bin=%d\n", n+1, regions->left_bin[n], regions->right_bin[n]);
    }
    if(regions->frac_defined[n]) {
      for(i = 0; i < verbose.indent; i++)
	printf(" ");
      printf("  Region %d: left frac=%f, right frac=%f\n", n+1, regions->left_frac[n], regions->right_frac[n]);
    }
  }
}

/* Counts the uniqe bins specified by regions. It does not double
   count bins which appear in multiple selections. */
int count_nrbins_in_Regions(pulselongitude_regions_definition *regions, int debug)
{
  int i, j, k, total, ok;
  if(regions == NULL)
    return 0;
  total = 0;
  for(i = 0; i < regions->nrRegions; i++) {
    if(regions->bins_defined[i] == 0) {
      fflush(stdout);
      printwarning(debug, "WARNING: count_nrbins_in_Regions called without region defined in bins");
      return 0;
    }
    if(i == 0) {
      total += regions->right_bin[i] - regions->left_bin[i] + 1;
    }else {
      for(j = regions->left_bin[i]; j <= regions->right_bin[i]; j++) {
	ok = 1;
	for(k = 0; k < i; k++) {
	  if(j >= regions->left_bin[k] && j <= regions->right_bin[k])
	    ok = 0;
	}
	if(ok)
	  total ++;
      }
    }
  }
  return total;
}

//START REGION DEVELOP
//START REGION RELEASE

// Return 1 = success, 0 = fail
int initPulselongitudeRegion(pulselongitude_regions_definition *region, verbose_definition verbose)
{
  region->bins_defined = malloc(MAX_pulselongitude_regions*sizeof(int));
  region->left_bin = malloc(MAX_pulselongitude_regions*sizeof(int));
  region->right_bin = malloc(MAX_pulselongitude_regions*sizeof(int));
  region->frac_defined = malloc(MAX_pulselongitude_regions*sizeof(int));
  region->left_frac = malloc(MAX_pulselongitude_regions*sizeof(float));
  region->right_frac = malloc(MAX_pulselongitude_regions*sizeof(float));
  if(region->bins_defined == NULL || region->left_bin == NULL || region->right_bin == NULL || region->frac_defined == NULL || region->left_frac == NULL || region->right_frac == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR initPulselongitudeRegion: Memory allocation error");
    return 0;
  }
  clearPulselongitudeRegion(region);
  return 1;
}

// Release the memory associated with an onpulse region
void freePulselongitudeRegion(pulselongitude_regions_definition *region)
{
  free(region->bins_defined);
  free(region->left_bin);
  free(region->right_bin);
  free(region->frac_defined);
  free(region->left_frac);
  free(region->right_frac);
}

// Copy a region definition to an already initialised other onpulse region
void copyPulselongitudeRegion(pulselongitude_regions_definition source, pulselongitude_regions_definition *destination)
{
  int i;
  destination->nrRegions = source.nrRegions;
  if(source.nrRegions > 0) {
    for(i = 0; i < source.nrRegions; i++) {
      destination->bins_defined[i] = source.bins_defined[i];
      destination->left_bin[i] = source.left_bin[i];
      destination->right_bin[i] = source.right_bin[i];
      destination->frac_defined[i] = source.frac_defined[i];
      destination->left_frac[i] = source.left_frac[i];
      destination->right_frac[i] = source.right_frac[i];
    }
  }
}

// This function is automatically called from within initPulselongitudeRegion
void clearPulselongitudeRegion(pulselongitude_regions_definition *region)
{
  int i;
  region->nrRegions = 0;
  for(i = 0; i < MAX_pulselongitude_regions; i++) {
    region->frac_defined[i] = 0;
    region->bins_defined[i] = 0;
  }
}

//START REGION DEVELOP

/* Given the array selection of nrbins integers, construct regions
   based on those bins where selection != 0. Returns 1 on success, 0
   on error. */
int selection_to_region(int *selection, int nrbins, pulselongitude_regions_definition *onpulse, verbose_definition verbose)
{
  int i, status;
  clearPulselongitudeRegion(onpulse);
  status = 0;
  for(i = 0; i < nrbins; i++) {
    if(selection[i] && status == 0) {  // Start new region
      if(onpulse->nrRegions == MAX_pulselongitude_regions) {
	fflush(stdout);
	printerror(verbose.debug, "ERROR selection_to_region: Max nr of regions exceeded");
	//	int j;
	//	for(j = 0; j < nrbins; j++) {
	//	  printerror(verbose.debug, "DEBUG: %d %d", j, selection[j]);
	//	}
	return 0;
      }
      onpulse->bins_defined[onpulse->nrRegions] = 1;
      onpulse->left_bin[onpulse->nrRegions] = i;
      status = 1;
    }else if(selection[i] && status == 1) {  // Continuing new region
    }else if(selection[i] == 0 && status == 1) {  // Closing new region
      onpulse->right_bin[onpulse->nrRegions] = i-1;
      onpulse->nrRegions++;
      status = 0;
    }else if(selection[i] == 0 && status == 0) {  // Continuing not selected region
    }else {
      fflush(stdout);
      printerror(verbose.debug, "ERROR selection_to_region: This shouldn't happen");
      return 0;
    }
  }
  if(status == 1) {                    // Close last opened region
    onpulse->right_bin[onpulse->nrRegions] = nrbins-1;
    onpulse->nrRegions++;
    status = 0;
  }
  return 1;
}

//START REGION DEVELOP
//START REGION RELEASE


/* Convert all regions defined as fracions to bin numbers via: 
   int = frac*scale + offset
*/
void region_frac_to_int(pulselongitude_regions_definition *region, float scale, float offset)
{
  int i;
  for(i = 0; i < region->nrRegions; i++) {
    if(region->frac_defined[i]) {
      region->left_bin[i] = region->left_frac[i]*scale + offset;
      region->right_bin[i] = region->right_frac[i]*scale + offset;
      region->bins_defined[i] = 1;
    }
  }
}

//START REGION DEVELOP
//START REGION RELEASE

/* Convert all regions defined as bin numbers to pulse phase via: 
   frac = bin*scale + offset
*/
void region_int_to_frac(pulselongitude_regions_definition *region, float scale, float offset)
{
  int i;
  for(i = 0; i < region->nrRegions; i++) {
    if(region->bins_defined[i]) {
      region->frac_defined[i] = 1;
      region->left_frac[i] = region->left_bin[i]*scale + offset;
      region->right_frac[i] = region->right_bin[i]*scale + offset;
    }
  }
}

//START REGION DEVELOP
//START REGION RELEASE

/* Show the selected regions as command line options (option). Where
   could be set to stdout for instance. If optionFrac != NULL and
   region contains a definition in terms of pulse phase, that will be
   outputted as well. */
void regionShowNextTimeUse(pulselongitude_regions_definition region, char *option, char *optionFrac, FILE *where)
{
  int i, ok;
  if(region.nrRegions > 0) {
    fprintf(where, "If you repeat this command, you could specify on the command line: ");
    for(i = 0; i < region.nrRegions; i++) {
      if(region.bins_defined[i])
	fprintf(where, "%s '%d %d' ", option, region.left_bin[i], region.right_bin[i]);
    }
    if(optionFrac != NULL) {
      ok = 0;
      for(i = 0; i < region.nrRegions; i++) {
	if(region.frac_defined[i])
	  ok = 1;
      }
      if(ok) {
	fprintf(where, "\n  alternatively you can use the command line:    ");
	for(i = 0; i < region.nrRegions; i++) {
	  if(region.frac_defined[i])
	    fprintf(where, "%s '%f %f' ", optionFrac, region.left_frac[i], region.right_frac[i]);
	  else if(region.bins_defined[i])
	    fprintf(where, "%s '%d %d' ", option, region.left_bin[i], region.right_bin[i]);
	}
      }
    }
    fprintf(where, "\n");
  }
}

//START REGION DEVELOP


int getstrings(char *txt, int nrstrings, char *string1, char *string2, char *string3, int debug)
{
  int ret;
  char *txtptr, *txtptr2, *txtptr3;
  ret = 0;
  txtptr = txt;

  if(nrstrings >= 1) {
    txtptr2 = strstr(txtptr, "'");
    if(txtptr2 == NULL) {
      fflush(stdout);
      printerror(debug, "ERROR getstrings: Cannot parse '%s'", txt);
      return ret;
    }
    txtptr3 = strstr(txtptr2+1, "'");
    if(txtptr3 == NULL) {
      fflush(stdout);
      printerror(debug, "ERROR getstrings: Cannot parse '%s'", txt);
      return ret;
    }
    /*    printf("Found substring \"%s\"\n", txtptr2);
    printf("Found substring \"%s\"\n", txtptr3);
    */
    strncpy(string1, txtptr2+1, txtptr3-txtptr2);
    string1[txtptr3-txtptr2-1] = 0;
    /*    printf("Found string \"%s\"\n", string1); */
    ret++;
    nrstrings--;
    txtptr = txtptr3+1;
  }
  if(nrstrings >= 1) {
    txtptr2 = strstr(txtptr, "'");
    if(txtptr2 == NULL) {
      fflush(stdout);
      printerror(debug, "ERROR getstrings: Cannot parse '%s'", txt);
      return ret;
    }
    txtptr3 = strstr(txtptr2+1, "'");
    if(txtptr3 == NULL) {
      fflush(stdout);
      printerror(debug, "ERROR getstrings: Cannot parse '%s'", txt);
      return ret;
    }
    /*    printf("Found substring \"%s\"\n", txtptr2);
    printf("Found substring \"%s\"\n", txtptr3);
    */
    strncpy(string2, txtptr2+1, txtptr3-txtptr2);
    string2[txtptr3-txtptr2-1] = 0;
    /*    printf("Found string \"%s\"\n", string2); */
    ret++;
    nrstrings--;
    txtptr = txtptr3+1;
  }
  if(nrstrings >= 1) {
    txtptr2 = strstr(txtptr, "'");
    if(txtptr2 == NULL) {
      fflush(stdout);
      printerror(debug, "Cannot parse '%s'", txt);
      return ret;
    }
    txtptr3 = strstr(txtptr2+1, "'");
    if(txtptr3 == NULL) {
      fflush(stdout);
      printerror(debug, "ERROR getstrings: Cannot parse '%s'", txt);
      return ret;
    }
    /*    printf("Found substring \"%s\"\n", txtptr2);
    printf("Found substring \"%s\"\n", txtptr3);
    */
    strncpy(string3, txtptr2+1, txtptr3-txtptr2);
    string3[txtptr3-txtptr2-1] = 0;
    /*    printf("Found string \"%s\"\n", string3); */
    ret++;
    nrstrings--;
    txtptr = txtptr3+1;
  }

  return ret;
}

int getfloat(char *txt, int c, float *f)
{
  long i, j;
  char txt2[PgplotLibMaxStringLength], *txt3;
  strcpy(txt2, txt);
  // Replace tabs by spaces
  for(i = 0; i < strlen(txt2); i++) {
    if(txt2[i] == '\t')
      txt2[i] = ' ';
  }
  i = 0;
  do {
    if(i == 0) {
      txt3 = strtok(txt2, " ");
    }else {
      txt3 = strtok(NULL, " ");
    }
    if(txt3 != NULL) {
      i++;
      if(i == c) {
	j = sscanf(txt3, "%f", f);
	if(j != 1)
	  return 0;
      }
    }
  }while(txt3 != NULL && i < c);
  if(txt3 == NULL)
    return 0;
  return 1;
}

int getdouble(char *txt, int c, double *f)
{
  long i, j;
  char txt2[PgplotLibMaxStringLength], *txt3;
  strcpy(txt2, txt);
  // Replace tabs by spaces
  for(i = 0; i < strlen(txt2); i++) {
    if(txt2[i] == '\t')
      txt2[i] = ' ';
  }
  i = 0;
  do {
    if(i == 0) {
      txt3 = strtok(txt2, " ");
    }else {
      txt3 = strtok(NULL, " ");
    }
    if(txt3 != NULL) {
      i++;
      if(i == c) {
	j = sscanf(txt3, "%lf", f);
	if(j != 1)
	  return 0;
      }
    }
  }while(txt3 != NULL && i < c);
  if(txt3 == NULL)
    return 0;
  return 1;
}

int getfloats(FILE *fin, float *f1, int c1, float *f2, int c2, float *f3, int c3)
{
  int j;
  char txt[PgplotLibMaxStringLength];
  if(fgets(txt, PgplotLibMaxStringLength, fin) == NULL)
    return 0;
  else {
    j = 0;
    if(c1 > 0)
      j += getfloat(txt, c1, f1);
    if(c2 > 0)
      j += getfloat(txt, c2, f2);
    if(c3 > 0)
      j += getfloat(txt, c3, f3);
  }
  return j;
}


void printPgplotInterpratorHelp()
{
  printf("Supported device commands\n");
  printf("  pgopen    devicename\n");
  printf("  pgbbuf\n");
  printf("  pgebuf\n");
  printf("  pgclos\n");
  printf("  pgend\n");
  printf("  pgpap     width[float] aspect[float]\n");
  printf("  pgask     flag[int]\n");
  printf("  pgpage  \n");
  printf("  size      x[int] y[int]\n");
  printf("\nSupported viewport commands\n");
  printf("  pgsvp     left[float] right[float] bottom[float] top[float]\n");
  printf("  pgswin    xmin[float] xmax[float] ymin[float] ymax[float]\n");
  printf("  pglab     'XLBL[sting]' 'YLBL[sting]' 'TOPLBL[sting]'\n");
  printf("  pgbox     xopt[sting] xtick[float] nxsub[int] yopt[sting] ytick[float] nysub[int]\n");
  printf("  pgaxis    'opt[string]' x1[float] y1[float] x2[float] y2[float] v1[float] v2[float] step[float]\n");
  printf("            nsub[int] dmajl[float] dmajr[float] fmin[float] disp[float] orient[float]\n");
  printf("  pgtick    x1[float] y1[float] x2[float] y2[float] v[float] tikl[float] rikr[float] disp[float]\n");
  printf("            orient[float] 'str[string]'\n");
  printf("  pgwedg    side[string] disp[float] width[float] fg[float] bg[float] 'label'\n");
  printf("  pgctab    filename nc[int] contra[float] bright[float]\n");
  printf("            binary file should contain nc L values, then nc R values, nc G values and nc B values.\n");
  printf("\nSupported pen commands\n");
  printf("  pgctab    file[string] nc[int] contra[float] bright[float]\n");
  printf("  pgscf     font[int] (normal/roman/italic/script)\n");
  printf("  pgsch     characterheight[float]\n");
  printf("  pgsci     colorindex[int]\n");
  printf("  pgscir    icilo[int] icihi[int]\n");
  printf("  pgscr     colorindex[int] cr[float] cg[float] cb[float]\n");
  printf("  pgsfs     fillstyle[int] (solid,outline,hatched,cross-hatched).\n");
  printf("  pgshs     angle[float] sepn [float] phase[float]\n");
  printf("  pgsitf    itf[int]\n");
  printf("  pgsls     lynestyle[int] (full,dashed,dot-dash,dotted,dash-dot-dot-dot).\n"); 
  printf("  pgslw     linewidth[int]\n");
  printf("\nSupported text commands\n");
  printf("  pgmtxt      side[string] disp[float] coord[float] fjust[float] 'text[string]'   (text along axis)\n");
  printf("  pgptxt      x[float] y[float] angle[float] fjust[float] 'text[string]'\n");
  printf("  pgtext      x[float] y[float] 'text[string]'\n");
  printf("  textoutline colorinside[int] coloroutside[int] linewidthinside[int] linewidthrim[int]\n");
  printf("              Set linewidthrim=0 to disable\n");
  printf("\nSupported draw commands\n");
  printf("  pgarro    x1[float] y1[float] x2[float] y2[float]\n");
  printf("  pgcirc    xcent[float] ycent[float] radius[float]\n");
  printf("  pgdraw    x[float] y[float]\n");
  printf("  pgerr1    dir[int] x[float] y[float] err[float] t[float]\n");
  printf("  pggray    file[string] idim[int] jdim[int] i1[int] i2[int] j1[int] j2[int] fg[float] bg[float]\n");
  printf("            tr0[float] tr1[float] tr2[float] tr3[float] tr4[float] tr5[float]\n");
  printf("  pgimag    file[string] idim[int] jdim[int] i1[int] i2[int] j1[int] j2[int] a1[float] a2[float]\n");
  printf("            tr0[float] tr1[float] tr2[float] tr3[float] tr4[float] tr5[float]\n");
  printf("  pgcont    file[string] idim[int] jdim[int] i1[int] i2[int] j1[int] j2[int] nc[int] c[nc floats]\n");
  printf("            tr0[float] tr1[float] tr2[float] tr3[float] tr4[float] tr5[float]\n");
  printf("  pgconl    file[string] idim[int] jdim[int] i1[int] i2[int] j1[int] j2[int] c[floats] tr0[float]\n");
  printf("            tr1[float] tr2[float] tr3[float] tr4[float] tr5[float] intval[int] minint[int] 'label'\n");
  printf("  pgmove    x[float] y[float]\n");
  printf("  pgpoly    nrpoints[int] filename (binary floats)\n");
  printf("  pgpt1     x[float] y[float] pointtype[int]\n");
  printf("  point     x[float] y[float]       (use ptype to set point type)\n");
  printf("  pgrect    x1[float] y1[float] x2[float] y2[float]\n");
  printf("  area      color fillstyle nrpoints x1 y1 x2 y2 ...\n");
  printf("  line      x1 y1 x2 y2 linewidth linestyle\n");
  printf("  load      arrow       file[string] x1:y1:x2:y2\n");
  printf("            curve       file[string] x:y\n");
  printf("            hist        file[string] x:y\n");
  printf("            histfilled  file[string] x:y\n");
  printf("            xerrorbar   file[string] x:y:dx\n");
  printf("            xerrorbar2  file[string] x:y:xleft:xright\n");
  printf("            yerrorbar   file[string] x:y:dy\n");
  printf("            yerrorbar2  file[string] x:y:ybottom:ytop\n");
  printf("            point       file[string] x:y    (use ptype to set point type)\n");
  printf("            map         file[string] nx[int] ny[int]\n");
  printf("            mapA          same, but now file is an ascii file with numbers rather than binary map\n");
  printf("            mapAScale   file[string] nx[int] ny[int] TR0 TR1 TR2 TR3 TR4 TR5\n");
  printf("  ptype     pointtype[int] OR\n");
  printf("            pointtype1[int] color1[int] pointtype2[int] color2[int] (Can do outlined symbols)\n");
  printf("  var       %%number value  (to be used in load command)\n");
  printf("\nOther commands\n");
  printf("  help      this help, can also give a function to get more info\n");
  printf("  include   include this filename with pgplot commands\n");
  printf("  quit      stop reading from input\n");
}


void doPgplotInterprator_drawpoint(float x, float y, int pointtype, int pointtype_double, int pointtype2, int pointtype_c, int pointtype_c2)
{
  if(pointtype_double == 0) {
    ppgpt1(x, y, pointtype);
  }else {
    ppgsci(pointtype_c);
    ppgpt1(x, y, pointtype);
    ppgsci(pointtype_c2);
    ppgpt1(x, y, pointtype2);
  }
}


/* Executes pgplotinterprator command.
Return values:
 0   = error
 1   = processed ok
 100 = Reached quit command
*/
int doPgplotInterprator(char *txt, int debug)
{
  static int firsttime = 0, pointtype, pointtype2, pointtype_c, pointtype_c2, pointtype_double;
  static int nrvariables2 = 0;
  /*
  static float dx, dy, dx2, dy2, dx3, dy3;
  static float sx, sy, sx2, sy2, sx3, sy3;
  */
  int i, j, k, i1, i2, i3, i4, i5, i6, i7, i8, pointnr, nrplotvariables, continuerunning;
  static int textoutline_c1, textoutline_c2, textoutline_lw1, textoutline_lw2;
  /*, loadtype, nrcols */
  char pgcommand[PgplotLibMaxStringLength], plotvariables[5][PgplotLibMaxStringLength];
  char txt3[PgplotLibMaxStringLength], string1[PgplotLibMaxStringLength], string2[PgplotLibMaxStringLength], string3[PgplotLibMaxStringLength];
  char *tmpcharptr;
  float f[100];
  float f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, *map, *datax, *datay, tr[6], max, min;
  float x, y, x2, y2, oldf1, oldf2, olddx;  //, oldf3
  float xarr[100], yarr[100];
  static char variables[100][100];
  static double variables_values[100];
  double answer;
  int nrvariables;
  long ret;
  FILE *fin2;
  if(firsttime == 0) {
    printf("doPgplotInterprator: initializing variables\n");
    firsttime = 1;
    /*
    dx = dy = dx2 = dy2 = dx3 = dy3 = 0;
    sx = sy = sx2 = sy2 = sx3 = sy3 = 1;
    */
    pointtype = 1;
    pointtype_double = 0;
    textoutline_c1 = 1;
    textoutline_c2 = 1;
    textoutline_lw1 = 1;
    textoutline_lw2 = 0;
  }
  j = sscanf(txt, "%s", pgcommand);
  if(j == 1) {
    if(strcmp(pgcommand, "help") == 0) {
      j = sscanf(txt, "%s %s", pgcommand, string1);
      if(j == 2) {
	printSubjectHelp(string1);
      }else {
	printPgplotInterpratorHelp();
      }
    }else if(strcmp(pgcommand, "include") == 0) {
      sscanf(txt, "%s %s", pgcommand, string1);
      printf("Opening file '%s'\n", string1);
      fin2 = fopen(string1, "r");
      if(fin2 == NULL) {
	fflush(stdout);
	printerror(debug, "ERROR doPgplotInterprator: Cannot open '%s'", string1);
	return 0;
      }else {
	continuerunning = 1;
	do {
	  if(fgets(string1, PgplotLibMaxStringLength, fin2) == NULL) {
	    continuerunning = 0;
	  }else {
	    if(doPgplotInterprator(string1, debug) == 0) {
	      return 0;
	    }
	  }
	}while(continuerunning);
	fclose(fin2);
      }
    }else if(strcmp(pgcommand, "quit") == 0) {
      return 100;
    }else if(strcmp(pgcommand, "pgopen") == 0) {
      sscanf(txt, "%s %s", pgcommand, string1);
      printf("Opening device '%s'\n", string1);
      ppgopen(string1);
    }else if(strcmp(pgcommand, "pgclos") == 0) {
      ppgclos();
    }else if(strcmp(pgcommand, "pgend") == 0) {
      ppgend();
    }else if(strcmp(pgcommand, "pgbbuf") == 0) {
      ppgbbuf();
    }else if(strcmp(pgcommand, "pgebuf") == 0) {
      ppgebuf();
    }else if(strcmp(pgcommand, "pgsvp") == 0) {
      sscanf(txt, "%s %f %f %f %f", pgcommand, &f1, &f2, &f3, &f4);
      ppgsvp(f1, f2, f3, f4);
      /*      xleftvp = f1; 
      xrightvp = f2; 
      ybottomvp = f3;
      ytopvp = f4; */
    }else if(strcmp(pgcommand, "pgswin") == 0) {
      sscanf(txt, "%s %f %f %f %f", pgcommand, &f1, &f2, &f3, &f4);
      ppgswin(f1, f2, f3, f4);
      /*      xleft = f1; 
      xright = f2; 
      ybottom = f3;
      ytop = f4; */
    }else if(strcmp(pgcommand, "var") == 0) {
      sscanf(txt, "%s %s %lf", pgcommand, string1, &answer);
      variables_values[nrvariables2] = answer;
      strcpy(variables[nrvariables2++], string1);
      printf("Adding variable %s (%f)\n", string1, answer);
      if(string1[0] != '%') {
	fflush(stdout);
	printerror(debug, "ERROR doPgplotInterprator: Variable name should start with a percent sign instead of '%c'", string1[0]);
      }
    }else if(strcmp(pgcommand, "pgbox") == 0) {
      sscanf(txt, "%s %s %f %d %s %f %d", pgcommand, string1, &f1, &i1, string2, &f2, &i2);
      ppgbox(string1, f1, i1, string2, f2, i2);
    }else if(strcmp(pgcommand, "pgaxis") == 0) {
      sscanf(txt, "%s %s %f %f %f %f %f %f %f %d %f %f %f %f %f", pgcommand, string1, &f1, &f2, &f3, &f4, &f5, &f6, &f7, &i1, &f8, &f9, &f10, &f11, &f12);
      ppgaxis(string1, f1, f2, f3, f4, f5, f6, f7, i1, f8, f9, f10, f11, f12);
    }else if(strcmp(pgcommand, "pgtick") == 0) {
      sscanf(txt, "%s %f %f %f %f %f %f %f %f %f", pgcommand, &f1, &f2, &f3, &f4, &f5, &f6, &f7, &f8, &f9);
      if(getstrings(txt, 1, string1, string2, string3, debug) != 1) {
	fflush(stdout);
	printerror(debug, "ERROR doPgplotInterprator: Cannot parse '%s'", txt);
      }
      ppgtick(f1, f2, f3, f4, f5, f6, f7, f8, f9, string1);
    }else if(strcmp(pgcommand, "pgwedg") == 0) {
      sscanf(txt, "%s %s %f %f %f %f", pgcommand, string1, &f1, &f2, &f3, &f4);
      if(getstrings(txt, 1, string2, string3, string3, debug) != 1) {
	fflush(stdout);
	printerror(debug, "ERROR doPgplotInterprator: Cannot parse '%s'", txt);
      }
      ppgwedg(string1, f1, f2, f3, f4, string2);
    }else if(strcmp(pgcommand, "pglab") == 0) {
      if(getstrings(txt, 3, string1, string2, string3, debug) != 3) {
	fflush(stdout);
	printerror(debug, "ERROR doPgplotInterprator: Cannot parse '%s'", txt);
      }else
	ppglab(string1, string2, string3);
    }else if(strcmp(pgcommand, "pgpoly") == 0) {
      sscanf(txt, "%s %d %s", pgcommand, &i1, string1);
      datax = (float *)malloc(sizeof(float)*i1);
      datay = (float *)malloc(sizeof(float)*i1);
      if(datax == NULL || datay == NULL) {
	fflush(stdout);
	printerror(debug, "ERROR doPgplotInterprator: Cannot allocate memory");
	return 0;
      }
      fin2 = fopen(string1, "rb");
      if(fin2 == NULL) {
	fflush(stdout);
	printerror(debug, "ERROR doPgplotInterprator: Cannot open %s", string1);
	return 0;
      }
      ret = fread(datax, sizeof(float), i1, fin2);
      if(ret != i1) {
	fflush(stdout);
	printerror(debug, "ERROR doPgplotInterprator: Read error %s (%ld != %ld)", string1, ret, (long)i1);
	return 0;
      }
      ret = fread(datay, sizeof(float), i1, fin2);
      if(ret != i1) {
	fflush(stdout);
	printerror(debug, "ERROR doPgplotInterprator: Read error %s (%ld != %ld)", string1, ret, (long)i1);
	return 0;
      }
      ppgpoly(i1, datax, datay);
      free(datax);
      free(datay);
      fclose(fin2);
    }else if(strcmp(pgcommand, "pggray") == 0 || strcmp(pgcommand, "pgimag") == 0) {
      sscanf(txt, "%s %s %d %d %d %d %d %d %f %f %f %f %f %f %f %f", pgcommand, string1, &i1, &i2, &i3, &i4, &i5, &i6, &f1, &f2, &f[0], &f[1], &f[2], &f[3], &f[4], &f[5]);
      map = (float *)malloc(sizeof(float)*i1*i2);
      if(map == NULL) {
	fflush(stdout);
	printerror(debug, "ERROR doPgplotInterprator: Cannot allocate memory");
	return 0;
      }
      fin2 = fopen(string1, "rb");
      if(fin2 == NULL) {
	fflush(stdout);
	printerror(debug, "ERROR doPgplotInterprator: Cannot open %s", string1);
	return 0;
      }
      ret = fread(map, sizeof(float), i1*i2, fin2);
      if(ret != i1*i2) {
	fflush(stdout);
	printerror(debug, "ERROR doPgplotInterprator: Read error %s (%ld != %ld)", string1, ret, (long)i1*i2);
	return 0;
      }
      /*
      for(ret = 1; ret < i1*i2; ret++) {
	if(map[ret] != map[ret-1])
	  printf("%e %d %d\n", map[ret], i1, i2);
	  }*/
      if(strcmp(pgcommand, "pggray") == 0)
	ppggray(map, i1, i2, i3, i4, i5, i6, f1, f2, f);
      else
	ppgimag(map, i1, i2, i3, i4, i5, i6, f1, f2, f);
      /*  printf("  pggray  file[string] idim[int] jdim[int] i1[int] int[int] j1[int] j2[int] fg[float] bg[float] tr0[float] tr1[float] tr2[float] tr3[float] tr4[float] tr5[float])\n");*/
      free(map);
      fclose(fin2);
    }else if(strcmp(pgcommand, "pgctab") == 0) {
      sscanf(txt, "%s %s %d %f %f", pgcommand, string1, &i1, &f1, &f2);
      map = (float *)malloc(sizeof(float)*i1*4);
      if(map == NULL) {
	fflush(stdout);
	printerror(debug, "ERROR doPgplotInterprator: Cannot allocate memory");
	return 0;
      }
      fin2 = fopen(string1, "rb");
      if(fin2 == NULL) {
	fflush(stdout);
	printerror(debug, "ERROR doPgplotInterprator: Cannot open %s", string1);
	return 0;
      }
      ret = fread(map, sizeof(float), i1*4, fin2);
      if(ret != i1*4) {
	fflush(stdout);
	printerror(debug, "ERROR doPgplotInterprator: Read error %s (%ld != %ld)", string1, ret, (long)i1*4);
	return 0;
      }
      ppgctab(&map[0], &map[1*i1], &map[2*i1], &map[3*i1], i1, f1, f2);
      free(map);
      fclose(fin2);
    }else if(strcmp(pgcommand, "pgcont") == 0) {
      sscanf(txt, "%s %s %d %d %d %d %d %d %d", pgcommand, string1, &i1, &i2, &i3, &i4, &i5, &i6, &i7);
      strcpy(string2, txt);
      tmpcharptr = strtok(string2, " ");
      if(tmpcharptr == NULL) {
	fflush(stdout);
	printerror(debug, "ERROR doPgplotInterprator: Parse error in pgcont (1)");
	return 0;
      }
      for(i = 0; i < 8; i++) {
	tmpcharptr = strtok(NULL, " ");
	if(tmpcharptr == NULL) {
	  fflush(stdout);
	  printerror(debug, "ERROR doPgplotInterprator: Parse error in pgcont (2)");
	  return 0;
	}
      }
      for(i = 0; i < abs(i7); i++) {
	tmpcharptr = strtok(NULL, " ");
	if(tmpcharptr == NULL) {
	  fflush(stdout);
	  printerror(debug, "ERROR doPgplotInterprator: Parse error in pgcont (3)");
	  return 0;
	}
	j = sscanf(tmpcharptr, "%f", &f[6+i]);
	if(j != 1) {
	  fflush(stdout);
	  printerror(debug, "ERROR doPgplotInterprator: Parse error in pgcont");
	  return 0;
	}
	/*		  printf("i=%f\n", f[6+i]); */
      }
      for(i = 0; i < 6; i++) {
	tmpcharptr = strtok(NULL, " ");
	if(tmpcharptr == NULL) {
	  fflush(stdout);
	  printerror(debug, "ERROR doPgplotInterprator: Parse error in pgcont (4)");
	  return 0;
	}
	j = sscanf(tmpcharptr, "%f", &f[i]);
	if(j != 1) {
	  fflush(stdout);
	  printerror(debug, "ERROR doPgplotInterprator: Parse error in pgcont");
	  return 0;
	}
      }
      map = (float *)malloc(sizeof(float)*i1*i2);
      if(map == NULL) {
	fflush(stdout);
	printerror(debug, "ERROR doPgplotInterprator: Cannot allocate memory");
	return 0;
      }
      fin2 = fopen(string1, "rb");
      if(fin2 == NULL) {
	fflush(stdout);
	printerror(debug, "ERROR doPgplotInterprator: Cannot open %s", string1);
	return 0;
      }
      ret = fread(map, sizeof(float), i1*i2, fin2);
      if(ret != i1*i2) {
	fflush(stdout);
	printerror(debug, "ERROR doPgplotInterprator: Read error %s (%ld != %ld)", string1, ret, (long)i1*i2);
	return 0;
      }
      ppgcont(map, i1, i2, i3, i4, i5, i6, &f[6], i7, f);
      free(map);
      fclose(fin2);
    }else if(strcmp(pgcommand, "pgconl") == 0) {
      sscanf(txt, "%s %s %d %d %d %d %d %d %f %f %f %f %f %f %f %d %d", pgcommand, string1, &i1, &i2, &i3, &i4, &i5, &i6, &f1, &f[0], &f[1], &f[2], &f[3], &f[4], &f[5], &i7, &i8);
      map = (float *)malloc(sizeof(float)*i1*i2);
      if(map == NULL) {
	fflush(stdout);
	printerror(debug, "ERROR doPgplotInterprator: Cannot allocate memory");
	return 0;
      }
      fin2 = fopen(string1, "rb");
      if(fin2 == NULL) {
	fflush(stdout);
	printerror(debug, "ERROR doPgplotInterprator: Cannot open %s", string1);
	return 0;
      }
      ret = fread(map, sizeof(float), i1*i2, fin2);
      if(ret != i1*i2) {
	fflush(stdout);
	printerror(debug, "ERROR doPgplotInterprator: Read error %s (%ld != %ld)", string1, ret, (long)i1*i2);
	return 0;
      }
      if(getstrings(txt, 1, string1, string2, string3, debug) != 1) {
	fflush(stdout);
	printerror(debug, "Cannot parse '%s'", txt);
      }
      ppgconl(map, i1, i2, i3, i4, i5, i6, f1, f, string1, i7, i8);
      free(map);
      fclose(fin2);
    }else if(strcmp(pgcommand, "pgpage") == 0) {
      ppgpage();
    }else if(strcmp(pgcommand, "pgask") == 0) {
      sscanf(txt, "%s %d", pgcommand, &i1);
      ppgask(i1);
    }else if(strcmp(pgcommand, "pgslw") == 0) {
      sscanf(txt, "%s %d", pgcommand, &i1);
      ppgslw(i1);
    }else if(strcmp(pgcommand, "pgsfs") == 0) {
      sscanf(txt, "%s %d", pgcommand, &i1);
      ppgsfs(i1);
    }else if(strcmp(pgcommand, "pgsls") == 0) {
      sscanf(txt, "%s %d", pgcommand, &i1);
      ppgsls(i1);
    }else if(strcmp(pgcommand, "pgsch") == 0) {
      sscanf(txt, "%s %f", pgcommand, &f1);
      ppgsch(f1);
    }else if(strcmp(pgcommand, "pgsci") == 0) {
      sscanf(txt, "%s %d", pgcommand, &i1);
      ppgsci(i1);
    }else if(strcmp(pgcommand, "pgsitf") == 0) {
      sscanf(txt, "%s %d", pgcommand, &i1);
      ppgsitf(i1);
    }else if(strcmp(pgcommand, "pgscir") == 0) {
      sscanf(txt, "%s %d %d", pgcommand, &i1, &i2);
      ppgscir(i1, i2);
    }else if(strcmp(pgcommand, "pgscf") == 0) {
      sscanf(txt, "%s %d", pgcommand, &i1);
      ppgscf(i1);
    }else if(strcmp(pgcommand, "pgscr") == 0) {
      sscanf(txt, "%s %d %f %f %f", pgcommand, &i1, &f1, &f2, &f3);
      ppgscr(i1, f1, f2, f3);
    }else if(strcmp(pgcommand, "pgshs") == 0) {
      sscanf(txt, "%s %f %f %f", pgcommand, &f1, &f2, &f3);
      ppgshs(f1, f2, f3);
    }else if(strcmp(pgcommand, "size") == 0) {
      sscanf(txt, "%s %d %d", pgcommand, &i1, &i2);
      f1 = i1*0.01175548589341692789994673739445429916373;
      f2 = (i2-1.0)/(float)i1;
      ppgpap(f1, f2);
    }else if(strcmp(pgcommand, "pgpap") == 0) {
      sscanf(txt, "%s %f %f", pgcommand, &f1, &f2);
      ppgpap(f1, f2);
    }else if(strcmp(pgcommand, "pgmove") == 0) {
      sscanf(txt, "%s %f %f", pgcommand, &f1, &f2);
      ppgmove(f1, f2);
    }else if(strcmp(pgcommand, "pgdraw") == 0) {
      sscanf(txt, "%s %f %f", pgcommand, &f1, &f2);
      ppgdraw(f1, f2);
    }else if(strcmp(pgcommand, "pgarro") == 0) {
      sscanf(txt, "%s %f %f %f %f", pgcommand, &f1, &f2, &f3, &f4);
      ppgarro(f1, f2, f3, f4);
    }else if(strcmp(pgcommand, "pgpt1") == 0) {
      sscanf(txt, "%s %f %f %d", pgcommand, &f1, &f2, &i1);
      ppgpt1(f1, f2, i1);
    }else if(strcmp(pgcommand, "point") == 0) {
      sscanf(txt, "%s %f %f", pgcommand, &f1, &f2);
      doPgplotInterprator_drawpoint(f1, f2, pointtype, pointtype_double, pointtype2, pointtype_c, pointtype_c2);
    }else if(strcmp(pgcommand, "pgerr1") == 0) {
      sscanf(txt, "%s %d %f %f %f %f", pgcommand, &i1, &f1, &f2, &f3, &f4);
      ppgerr1(i1, f1, f2, f3, f4);
    }else if(strcmp(pgcommand, "pgcirc") == 0) {
      sscanf(txt, "%s %f %f %f", pgcommand, &f1, &f2, &f3);
      ppgcirc(f1, f2, f3);
    }else if(strcmp(pgcommand, "ptype") == 0) {
      if(sscanf(txt, "%s %d %d %d %d", pgcommand, &i1, &i2, &i3, &i4) == 5) {
	pointtype = i1;
	pointtype_c = i2;
	pointtype2 = i3;
	pointtype_c2 = i4;
	pointtype_double = 1;
      }else if(sscanf(txt, "%s %d", pgcommand, &i1) == 2) {
	pointtype = i1;
	pointtype_double = 0;
      }else {
	fflush(stdout);
	printerror(debug, "ERROR doPgplotInterprator: Cannot parse '%s'. Expected 4 or 1 integers", txt);	
      }
    }else if(strcmp(pgcommand, "textoutline") == 0) {
      if(sscanf(txt, "%s %d %d %d %d", pgcommand, &i1, &i2, &i3, &i4) == 5) {
	textoutline_c1 = i1;
	textoutline_c2 = i2;
	textoutline_lw1 = i3;
	textoutline_lw2 = i4;
      }else {
	fflush(stdout);
	printerror(debug, "ERROR doPgplotInterprator: Cannot parse '%s'. Expected 4 integers", txt);	
      }
    }else if(strcmp(pgcommand, "pgtext") == 0) {
      sscanf(txt, "%s %f %f", pgcommand, &f1, &f2);
      if(getstrings(txt, 1, string1, string2, string3, debug) != 1) {
	fflush(stdout);
	printerror(debug, "ERROR doPgplotInterprator: Cannot parse '%s'", txt);
      }else {
	if(textoutline_lw2 == 0) {
	  ppgtext(f1, f2, string1);
	}else {
	  ppgslw(textoutline_lw1+textoutline_lw2);
	  ppgsci(textoutline_c2);
	  ppgtext(f1, f2, string1);
	  ppgslw(textoutline_lw1);
	  ppgsci(textoutline_c1);
	  ppgtext(f1, f2, string1);
	}
      }
    }else if(strcmp(pgcommand, "pgmtxt") == 0) {
      sscanf(txt, "%s %s %f %f %f", pgcommand, txt3, &f1, &f2, &f3);
      if(getstrings(txt, 1, string1, string2, string3, debug) != 1) {
	fflush(stdout);
	printerror(debug, "ERROR doPgplotInterprator: Cannot parse '%s'", txt);
      }else {
	if(textoutline_lw2 == 0) {
	  ppgmtxt(txt3, f1, f2, f3, string1);
	}else {
	  ppgslw(textoutline_lw1+textoutline_lw2);
	  ppgsci(textoutline_c2);
	  ppgmtxt(txt3, f1, f2, f3, string1);
	  ppgslw(textoutline_lw1);
	  ppgsci(textoutline_c1);
	  ppgmtxt(txt3, f1, f2, f3, string1);
	}
      }
    }else if(strcmp(pgcommand, "pgptxt") == 0) {
      sscanf(txt, "%s %f %f %f %f", pgcommand, &f1, &f2, &f3, &f4);
      if(getstrings(txt, 1, string1, string2, string3, debug) != 1) {
	fflush(stdout);
	printerror(debug, "ERROR doPgplotInterprator: Cannot parse '%s'", txt);
      }else {
	if(textoutline_lw2 == 0) {
	  ppgptxt(f1, f2, f3, f4, string1);
	}else {
	  ppgslw(textoutline_lw1+textoutline_lw2);
	  ppgsci(textoutline_c2);
	  ppgptxt(f1, f2, f3, f4, string1);
	  ppgslw(textoutline_lw1);
	  ppgsci(textoutline_c1);
	  ppgptxt(f1, f2, f3, f4, string1);
	}
      }
    }else if(strcmp(pgcommand, "line") == 0) {
      sscanf(txt, "%s %f %f %f %f %f %d", pgcommand, &f1, &f2, &f3, &f4, &f5, &i1);
      ppgslw(f5);
      ppgsls(i1);
      ppgmove(f1, f2);
      ppgdraw(f3, f4);
      ppgslw(1);
      ppgsls(1);
    }else if(strcmp(pgcommand, "area") == 0) {
      j = sscanf(txt, "%s %d %d %d", pgcommand, &i1, &i2, &i3);
      if(j != 4) {
	fflush(stdout);
	printerror(debug, "ERROR doPgplotInterprator: Cannot parse %s option, need 3+2n values.", pgcommand);
      }else {
	if(i3 == 2) {
	  j = sscanf(txt, "%s %d %d %d %f %f %f %f", pgcommand, &i1, &i2, &i3, &xarr[0], &yarr[0], &xarr[1], &yarr[1]);
	  if(j != 4+2*i3) {
	    fflush(stdout);
	    printerror(debug, "ERROR doPgplotInterprator: Cannot parse %s option, need %d values.", pgcommand, 3+2*i3);
	  }
	}else if(i3 == 3) {
	  j = sscanf(txt, "%s %d %d %d %f %f %f %f %f %f", pgcommand, &i1, &i2, &i3, &xarr[0], &yarr[0], &xarr[1], &yarr[1], &xarr[2], &yarr[2]);
	  if(j != 4+2*i3) {
	    fflush(stdout);
	    printerror(debug, "ERROR doPgplotInterprator: Cannot parse %s option, need %d values.", pgcommand, 3+2*i3);
	  }
	}else if(i3 == 4) {
	  j = sscanf(txt, "%s %d %d %d %f %f %f %f %f %f %f %f", pgcommand, &i1, &i2, &i3, &xarr[0], &yarr[0], &xarr[1], &yarr[1], &xarr[2], &yarr[2], &xarr[3], &yarr[3]);
	  if(j != 4+2*i3) {
	    fflush(stdout);
	    printerror(debug, "ERROR doPgplotInterprator: Cannot parse %s option, need %d values.", pgcommand, 3+2*i3);
	  }
	}else {
	  fflush(stdout);
	  printerror(debug, "ERROR doPgplotInterprator: This number of points in 'area' command are not supported.");
	}
	ppgsfs(i2);
	if(i2 == 3) {
	  ppgshs(45, 3, 0);
	  ppgslw(3);
	}
	ppgsci(i1);
	/*	printf("poly: %d %f %f\n", i3, xarr[2], yarr[2]); */
	ppgpoly(i3, xarr, yarr);
	ppgsfs(2);
	ppgpoly(i3, xarr, yarr);
	ppgsci(1);
	ppgsfs(1);
	ppgslw(1);
      }
      /*    }else if(strcmp(pgcommand, "LOAD") == 0) {
      sscanf(txt, "%s %s", pgcommand, string1);
      i1 = 1;
      i2 = 2;
      i3 = -1;
      loadtype = 0;
      if(strcmp(string1, "curve") == 0) {
	sscanf(txt, "%s %s %s %d %d", pgcommand, string1, string2, &i1, &i2); 
	nrcols = 2;
	loadtype = 1;
      }else if(strcmp(string1, "yerrorbar") == 0) {
	sscanf(txt, "%s %s %s %d %d %d", pgcommand, string1, string2, &i1, &i2, &i3);
	nrcols = 3;
	loadtype = 2;
      }else if(strcmp(string1, "point") == 0) {
	sscanf(txt, "%s %s %s %d %d", pgcommand, string1, string2, &i1, &i2);
	nrcols = 2;
	loadtype = 3;
      }else if(strcmp(string1, "map") == 0) {
	sscanf(txt, "%s %s %s %d %d", pgcommand, string1, string2, &i1, &i2);
	nrcols = -1;
	loadtype = 4;
      }else if(strcmp(string1, "hist") == 0 || strcmp(string1, "histfilled") == 0) {
	sscanf(txt, "%s %s %s %d %d", pgcommand, string1, string2, &i1, &i2);
	nrcols = 2;
	if(strcmp(string1, "hist") == 0) 
	  loadtype = 5;
	else
	  loadtype = 6;
      }else {
      fflush(stdout);
	printerror(debug, "doPgplotInterprator: Undefined load type: '%s'", string1);
	return 0;
      }
      if(loadtype == 4) 
	fin2 = fopen(string2, "rb");
      else
	fin2 = fopen(string2, "r");
      if(fin2 == NULL) {
	sprintf(txt3, "doPgplotInterprator: Cannot open '%s' ", string2);
	perror(txt3);
	return 0;
      }
      i = 1;
      pointnr = 0;
      do {
	if(nrcols == -1) {
	  if(i == 1) {
	    map = (float *)malloc(sizeof(float)*i1*i2);
	    if(map == NULL) {
      fflush(stdout);
	      printerror(debug, "doPgplotInterprator: Memory allocation error.");
	      return 0;
	    }
	  }
	  j = fread(&map[(i-1)*i1], sizeof(float), i1, fin2);
	  if(j != i1) {
      fflush(stdout);
	    printerror(debug, "doPgplotInterprator: Read error (%d != %d).", j, i1);
	    return 0;
	  }
	  j = nrcols;
	  if(i == i2) {
	    tr[0] = 0;
	    tr[1] = 1;
	    tr[2] = 0;
	    tr[3] = 0;
	    tr[4] = 0;
	    tr[5] = 1;
	    max = min = map[0];
	    for(k = 0; k < i2; k++) {
	      for(j = 0; j < i1; j++) {
		if(map[k*i1+j] > max)
		  max = map[k*i1+j];
		if(map[k*i1+j] < min)
		  min = map[k*i1+j];
	      }
	    }
	    ppggray(map, i1, i2,  1, i1, 1, i2, max, min, tr);
	    free(map);
	    j = 0;
	  }
	  i++;
	}else {
	  j = getfloats(fin2, &f1, i1, &f2, i2, &f3, i3);
	  if(j == nrcols) {
	    if(loadtype == 1) {
	      if(i == 1)
		ppgmove(sx3*(sx2*(sx*(f1+dx)+dx2)+dx3), sy3*(sy2*(sy*(f2+dy)+dy2)+dy3));
	      else {
		ppgdraw(sx3*(sx2*(sx*(f1+dx)+dx2)+dx3), sy3*(sy2*(sy*(f2+dy)+dy2)+dy3));
	      }
	    }else if(loadtype == 2) {
	      ppgmove(sx3*(sx2*(sx*(f1+dx)+dx2)+dx3), sy3*(sy2*(sy*(f2+dy)+dy2)+dy3));
	      if(f3 > 0)
		ppgerr1(6, sx3*(sx2*(sx*(f1+dx)+dx2)+dx3), sy3*(sy2*(sy*(f2+dy)+dy2)+dy3), f3, 3);
	    }else if(loadtype == 3) {
	      ppgpt1(sx3*(sx2*(sx*(f1+dx)+dx2)+dx3), sy3*(sy2*(sy*(f2+dy)+dy2)+dy3), pointtype);
	    }else if((loadtype == 5 || loadtype == 6) && pointnr != 0) {
	      x = sx3*(sx2*(sx*(oldf1+dx)+dx2)+dx3);
	      y = sy3*(sy2*(sy*(oldf2+dy)+dy2)+dy3);
	      x2 = sx3*(sx2*(sx*(f1+dx)+dx2)+dx3);
	      y2 = sy3*(sy2*(sy*(f2+dy)+dy2)+dy3);
	      olddx = (x2-x);
	      if(loadtype == 5) {
		ppgmove(x-0.5*(x2-x), y);
		ppgdraw(x+0.5*(x2-x), y);
		ppgdraw(x+0.5*(x2-x), y2);
	      }else {
		ppgrect(x-0.5*(x2-x), x+0.5*(x2-x), y, -1e10);
	      }
	    }
	    i++;
	  }
	  oldf1 = f1;	
	  oldf2 = f2;
	  oldf3 = f3;
	}
	pointnr++;
      }while(j == nrcols);
      if(nrcols != -1 && loadtype == 5 && pointnr != 0) {
	x2 = sx3*(sx2*(sx*(f1+dx)+dx2)+dx3);
	y2 = sy3*(sy2*(sy*(f2+dy)+dy2)+dy3);
	ppgmove(x2-0.5*(olddx), y2);
	ppgdraw(x2+0.5*(olddx), y2);
      }
      fclose(fin2);
      printf("doPgplotInterprator: Loading '%s' with %d points\n", string2, i);*/
    }else if(strcmp(pgcommand, "load") == 0) {
      sscanf(txt, "%*s %s %s", string1, string2);
      /* Skip two words load & type */
      if(strcmp(string1, "map") != 0 && strcmp(string1, "mapA") != 0  && strcmp(string1, "mapAScale") != 0) { 
	i2 = -1;
	i3 = 0;
	for(i1 = 0; i1 < strlen(txt); i1++) {
	  if(txt[i1] == ' ') {
	    if(i2 < 0) {
	      i3++;
	      if(i3 == 3)
		break;
	    }
	    i2 = 1;
	  }else {
	    i2 = -1;
	  }
	}
	//	printf("Hoi: '%s'\n", txt+i1);
	tmpcharptr = strtok(txt+i1, ":");
	nrplotvariables = 0;
	do {
	  if(tmpcharptr != NULL) {
	    strcpy(plotvariables[nrplotvariables++], tmpcharptr);
	  }
	  tmpcharptr = strtok(NULL, ":");
	}while(tmpcharptr != NULL);
	/* Put a $ in of the column nr if there only is a number */
	for(i1 = 0; i1 < nrplotvariables; i1++) {
	  i2 = 0;
	  for(i3 = 0; i3 < strlen(plotvariables[i1]); i3++)
	    if(plotvariables[i1][i3] == '$') {
	      i2 = 1;
	      break;
	    }
	  if(i2 == 0) {
	    i3 = sscanf(plotvariables[i1], "%d", &i2);
	    if(i3 != 1) {
	      fflush(stdout);
	      printerror(debug, "ERROR doPgplotInterprator: Cannot parse '%s' as an integer.", plotvariables[i1]);
	      return 0;
	    }else {
	      sprintf(plotvariables[i1], "$%d", i2);
	    }
	  }
	  //	  printf("Hoi: '%s'\n", plotvariables[i1]); 
	}
      }

      /* open file */
      if(strcmp(string1, "map") == 0) 
	fin2 = fopen(string2, "rb");
      else
	fin2 = fopen(string2, "r");
      if(fin2 == NULL) {
	sprintf(txt3, "ERROR doPgplotInterprator: Cannot open '%s' ", string2);
	perror(txt3);
	return 0;
      }

      if(strcmp(string1, "curve") == 0) {
	if(nrplotvariables != 2) {
	  fflush(stdout);
	  printerror(debug, "ERROR doPgplotInterprator: Expecting two variables");
	  return 0;
	}
      }else if(strcmp(string1, "yerrorbar") == 0) {
	if(nrplotvariables != 3) {
	  fflush(stdout);
	  printerror(debug, "ERROR doPgplotInterprator: Expecting three variables");
	  return 0;
	}
      }else if(strcmp(string1, "yerrorbar2") == 0) {
	if(nrplotvariables != 4) {
	  fflush(stdout);
	  printerror(debug, "ERROR doPgplotInterprator: Expecting four variables");
	  return 0;
	}
      }else if(strcmp(string1, "xerrorbar") == 0) {
	if(nrplotvariables != 3) {
	  fflush(stdout);
	  printerror(debug, "ERROR doPgplotInterprator: Expecting three variables");
	  return 0;
	}
      }else if(strcmp(string1, "xerrorbar2") == 0) {
	if(nrplotvariables != 4) {
	  fflush(stdout);
	  printerror(debug, "ERROR doPgplotInterprator: Expecting four variables");
	  return 0;
	}
      }else if(strcmp(string1, "arrow") == 0) {
	if(nrplotvariables != 4) {
	  fflush(stdout);
	  printerror(debug, "ERROR doPgplotInterprator: Expecting four variables");
	  return 0;
	}
      }else if(strcmp(string1, "point") == 0) {
	if(nrplotvariables != 2) {
	  fflush(stdout);
	  printerror(debug, "ERROR doPgplotInterprator: Expecting two variables");
	  return 0;
	}
	/* Counting variables does not appear to work for maps */
	/*      }else if(strcmp(string1, "map") == 0) {     
	if(nrplotvariables != 2) {
      fflush(stdout);
	  printerror(debug, "ERROR doPgplotInterprator: Expecting two variables");
	  return 0;
	}
      }else if(strcmp(string1, "mapA") == 0) {
	if(nrplotvariables != 2) {
      fflush(stdout);
	  printerror(debug, "ERROR doPgplotInterprator: Expecting two variables (got %d)", nrplotvariables);
	  return 0;
	  }  */
      }else if(strcmp(string1, "map") == 0) {     
      }else if(strcmp(string1, "mapA") == 0) {
      }else if(strcmp(string1, "mapAScale") == 0) {
      }else if(strcmp(string1, "hist") == 0 || strcmp(string1, "histfilled") == 0) {
	if(nrplotvariables != 2) {
	  fflush(stdout);
	  printerror(debug, "ERROR doPgplotInterprator: Expecting two variables");
	  return 0;
	}
      }else {
	fflush(stdout);
	printerror(debug, "ERROR doPgplotInterprator: Undefined load type: '%s'", string1);
	return 0;
      }

      i = 1;
      pointnr = 0;
      do {
	if(strcmp(string1, "map") == 0 || strcmp(string1, "mapA") == 0 || strcmp(string1, "mapAScale") == 0) {
	  if(i == 1) {
	    if(strcmp(string1, "mapAScale") == 0) {
	      sscanf(txt, "%*s %*s %*s %d %d %f %f %f %f %f %f", &i1, &i2, &f1, &f2, &f3, &f4, &f5, &f6);
	    }else {
	      f1 = 0;
	      f2 = 1;
	      f3 = 0;
	      f4 = 0;
	      f5 = 0;
	      f6 = 1;
	      sscanf(txt, "%*s %*s %*s %d %d", &i1, &i2);
	    }
	    map = (float *)malloc(sizeof(float)*i1*i2);
	    if(map == NULL) {
	      fflush(stdout);
	      printerror(debug, "ERROR doPgplotInterprator: Memory allocation error.");
	      return 0;
	    }
	  }
	  if(strcmp(string1, "map") == 0) {
	    j = fread(&map[(i-1)*i1], sizeof(float), i1, fin2);
	    if(j != i1) {
	      fflush(stdout);
	      printerror(debug, "ERROR doPgplotInterprator: Read error (%d != %d).", j, i1);
	      return 0;
	    }
	  }else {
	    for(j = 0; j < i1; j++) {
	      fscanf(fin2, "%f", &map[(i-1)*i1+j]);
	    }
	  }
	  if(i == i2) {
	    tr[0] = f1;
	    tr[1] = f2;
	    tr[2] = f3;
	    tr[3] = f4;
	    tr[4] = f5;
	    tr[5] = f6;
	    max = min = map[0];
	    for(k = 0; k < i2; k++) {
	      for(j = 0; j < i1; j++) {
		if(map[k*i1+j] > max)
		  max = map[k*i1+j];
		if(map[k*i1+j] < min)
		  min = map[k*i1+j];
	      }
	    }
	    ppggray(map, i1, i2,  1, i1, 1, i2, max, min, tr);
	    free(map);
	    pointnr = -999;
	  }
	  i++;
	}else { /* No map */
	  nrvariables = nrvariables2;
	  //	  for(i1 = 0; i1 < nrplotvariables; i1++) {
	  //	    printf("XXXX Expression%d = %s\n", i1+1, plotvariables[i1]);
	  //	  }
	  for(i1 = 0; i1 < nrplotvariables; i1++) {
	    strcpy(string3, plotvariables[i1]);
	    if(string3[0] != '$')
	      tmpcharptr = strtok(string3, "$");
	    else
	      tmpcharptr = string3;
	    do {
	      if(tmpcharptr != NULL) {
		if(tmpcharptr[0] == '$')
		  tmpcharptr = strtok(string3, "$");
		else
		  tmpcharptr = strtok(NULL, "$");
	      }
	      if(tmpcharptr != NULL) {
		k = sscanf(tmpcharptr, "%d", &i2);
		if(k != 1) {
		  /* Not an error anymore, assume it is a defined variable */
		  /*		        fflush(stdout);
printerror(debug, "ERROR doPgplotInterprator: Parse error in '%s'.", tmpcharptr);
		    return 0; */
		}else {
		  sprintf(variables[nrvariables++], "$%d", i2);
		  /*		  printf("Variables: '%s'\n", variables[nrvariables-1]); */
		}
	      }
	    }while(tmpcharptr != NULL);
	  }
	  //	  for(i1 = 0; i1 < nrvariables; i1++) {
	  //	    printf("XXXX Variable%d = %s\n", i1+1, variables[i1]);
	  //	  }
	  /* Now we know what varables to read, read a line from file. */
	  do {
	    tmpcharptr = fgets(string3, PgplotLibMaxStringLength, fin2);
	    if(tmpcharptr == NULL) {
	      pointnr = -999;
	      i2 = 1;
	      break;
	    }else {
	      i2 = 0;
	      for(i1 = 0; i1 < strlen(string3); i1++) {
		if(string3[i1] != ' ' && string3[i1] != '\n' && string3[i1] != '\r' && string3[i1] != '\t') {
		  if(string3[i1] != '#')
		    i2 = 1;
		  break;
		}
	      }
	    }
	  }while(i2 == 0);
	  if(pointnr >= 0) {
	    int allColsReadInOk;
	    /*	    printf("Accepting line: %s\n", string3); */
	    allColsReadInOk = 1;
	    for(i1 = nrvariables2; i1 < nrvariables; i1++) {
	      sscanf(variables[i1], "$%d", &i2);
	      i3 = getdouble(string3, i2, &variables_values[i1]);
	      if(i3 != 1) {
		fflush(stdout);
		printerror(debug, "ERROR doPgplotInterprator: error reading column %d in line %d", i2, pointnr+1);
		//		printerror(debug, "ERROR doPgplotInterprator: line reads \"%s\"", string3);
		allColsReadInOk = 0;
		/*		return 0;*/
	      }
	      /*	      printf("%s = %f\n",variables[i1], variables_values[i1]); */
	    }
	    //	    for(i1 = 0; i1 < nrvariables; i1++) {
	    //	      printf("XXXX Values%d = %e\n", i1+1, variables_values[i1]);
	    //	    }
	    if(allColsReadInOk) {
	      for(i1 = 0; i1 < nrplotvariables; i1++) {
		//		printf("Expr='%s'\n", plotvariables[i1]);
		//		for(i = 0; i < nrvariables; i++)
		//		  printf("  variable%d '%s' %e\n", i+1, variables[i], variables_values[i]);
		i2 = calc_expression(plotvariables[i1], nrvariables, variables, variables_values, &answer, 0);
		//		printf("  answer = %e (ret=%d)\n", answer, i2);
		if(i2 == 0) {
		  return 0;
		}
		f[i1] = answer;
	      }
	      /* blaat */
	      if(strcmp(string1, "curve") == 0) {
		if(i == 1)
		  ppgmove(f[0], f[1]);
		else {
		  ppgdraw(f[0], f[1]);
		}
	      }else if(strcmp(string1, "xerrorbar") == 0) {
		ppgmove(f[0], f[1]);
		ppgerr1(5, f[0], f[1], f[2], 3);
	      }else if(strcmp(string1, "yerrorbar") == 0) {
		ppgmove(f[0], f[1]);
		ppgerr1(6, f[0], f[1], f[2], 3);
	      }else if(strcmp(string1, "yerrorbar2") == 0) {
		ppgmove(f[0], f[1]);
		//		ppgerr1(6, f[0], f[1], f[1]-f[2], 4);
		//		ppgerr1(6, f[0], f[1], f[3]-f[1], 2);
		ppgerr1(4, f[0], f[1], f[1]-f[2], 3);
		ppgerr1(2, f[0], f[1], f[3]-f[1], 3);
	      }else if(strcmp(string1, "xerrorbar2") == 0) {
		ppgmove(f[0], f[1]);
		ppgerr1(3, f[0], f[1], f[1]-f[2], 3);
		ppgerr1(1, f[0], f[1], f[3]-f[1], 3);
	      }else if(strcmp(string1, "arrow") == 0) {
		ppgarro(f[0], f[1], f[2], f[3]);
	      }else if(strcmp(string1, "point") == 0) {
		doPgplotInterprator_drawpoint(f[0], f[1], pointtype, pointtype_double, pointtype2, pointtype_c, pointtype_c2);
	      }else if((strcmp(string1, "hist") == 0 || strcmp(string1, "histfilled") == 0) && pointnr != 0) {
		x = oldf1;
		y = oldf2;
		x2 = f[0];
		y2 = f[1];
		olddx = (x2-x);
		if(strcmp(string1, "hist") == 0) {
		  ppgmove(x-0.5*(x2-x), y);
		  ppgdraw(x+0.5*(x2-x), y);
		  ppgdraw(x+0.5*(x2-x), y2);
		}else {
		  ppgrect(x-0.5*(x2-x), x+0.5*(x2-x), y, -1e10);
		}
	      }
	      i++;
	      oldf1 = f[0];	
	      oldf2 = f[1];
	      //	      oldf3 = f[2];
	      pointnr++;
	    }
	  }
	}
      }while(pointnr >= 0);
      if(strcmp(string1, "hist") == 0) {
	x2 = f[0];
	y2 = f[1];
	ppgmove(x2-0.5*(olddx), y2);
	ppgdraw(x2+0.5*(olddx), y2);
      }
      fclose(fin2);
      printf("doPgplotInterprator: Loading '%s' with %d points\n", string2, i-1);
    }else {
      if(pgcommand[0] != '#')  {
	fflush(stdout);
	printerror(debug, "ERROR doPgplotInterprator: Unknown command: '%s'", pgcommand);
      }
    }
  }  // End of if(j == 1)
  return 1;
}


int commandline_PgplotInterprator(int argc, char **argv, int debug)
{
  int i;
  for(i = 0; i < argc; i++) {
    if(strcmp(argv[i], "-pgplot") == 0) {
      if(doPgplotInterprator(argv[i+1], debug) == 0)
	return 0;
    }
  }
  return 1;
}

//START REGION RELEASE

/* Show the -cmap options */
void printCMAPCommandlineOptions(FILE *printdevice)
{
  fprintf(printdevice, "Valid options for the -cmap option are:\n");
  fprintf(printdevice, "  GRAYSCALE\n");
  fprintf(printdevice, "  INVERTED_GRAYSCALE\n");
  fprintf(printdevice, "  RED\n");
  fprintf(printdevice, "  INVERTED_RED\n");
  fprintf(printdevice, "  GREEN\n");
  fprintf(printdevice, "  INVERTED_GREEN\n");
  fprintf(printdevice, "  BLUE\n");
  fprintf(printdevice, "  INVERTED_BLUE\n");
  fprintf(printdevice, "  CYAN\n");
  fprintf(printdevice, "  INVERTED_CYAN\n");
  fprintf(printdevice, "  HEAT\n");
  fprintf(printdevice, "  INVERTED_HEAT\n");
  fprintf(printdevice, "  HEAT2\n");
  fprintf(printdevice, "  INVERTED_HEAT2\n");
  fprintf(printdevice, "  HEAT3\n");
  fprintf(printdevice, "  INVERTED_HEAT3\n");
  fprintf(printdevice, "  HEAT4\n");
  fprintf(printdevice, "  INVERTED_HEAT4\n");
  fprintf(printdevice, "  COLD\n");
  fprintf(printdevice, "  INVERTED_COLD\n");
  fprintf(printdevice, "  DIVREDBLUE\n");
  fprintf(printdevice, "  INVERTED_DIVREDBLUE (a.k.a the Dutch flag)\n");
//START REGION DEVELOP
  fprintf(printdevice, "  PLASMA\n");
  fprintf(printdevice, "  INVERTED_PLASMA\n");
  fprintf(printdevice, "  FOREST\n");
  fprintf(printdevice, "  INVERTED_FOREST\n");
  fprintf(printdevice, "  ALIEN_GLOW\n");               // NOTE THAT THIS IS ALMOST THE SAME AS THE HEAT2 MAP
  fprintf(printdevice, "  INVERTED_ALIEN_GLOW\n");
//START REGION RELEASE
  fprintf(printdevice, "  INFERNO\n");
}


/* Returns 0 if parse error, or else it returns the map type.  */
int cmap_parse_commandline(int argc, char **argv, int debug)
{
  int i, j, type;
  char identifier[100];

  type = 0;

  for(i = 1; i < argc - 1; i++) {
    if(strcmp(argv[i], "-cmap") == 0) {
      i++;
      j = sscanf(argv[i], "%s", identifier);
      if(j != 1) {
	fflush(stdout);
	printerror(debug, "ERROR cmap_parse_commandline:  Cannot parse -cmap option.");
	printCMAPCommandlineOptions(stderr);
	return 0;
      }else {
	if(strcasecmp(identifier,"GRAYSCALE") == 0) {
	  type = PPGPLOT_GRAYSCALE;
	}else if(strcasecmp(identifier,"INVERTED_GRAYSCALE") == 0) {
	  type = PPGPLOT_INVERTED_GRAYSCALE;

	}else if(strcasecmp(identifier,"RED") == 0) {
	  type = PPGPLOT_RED;
	}else if(strcasecmp(identifier,"GREEN") == 0) {
	  type = PPGPLOT_GREEN;
	}else if(strcasecmp(identifier,"BLUE") == 0) {
	  type = PPGPLOT_BLUE;
	}else if(strcasecmp(identifier,"CYAN") == 0) {
	  type = PPGPLOT_CYAN;
	}else if(strcasecmp(identifier,"INVERTED_RED") == 0) {
	  type = PPGPLOT_INVERTED_RED;
	}else if(strcasecmp(identifier,"INVERTED_GREEN") == 0) {
	  type = PPGPLOT_INVERTED_GREEN;
	}else if(strcasecmp(identifier,"INVERTED_BLUE") == 0) {
	  type = PPGPLOT_INVERTED_BLUE;
	}else if(strcasecmp(identifier,"INVERTED_CYAN") == 0) {
	  type = PPGPLOT_INVERTED_CYAN;
	}else if(strcasecmp(identifier,"HEAT") == 0) {
	  type = PPGPLOT_HEAT;
	}else if(strcasecmp(identifier,"INVERTED_HEAT") == 0) {
	  type = PPGPLOT_INVERTED_HEAT;
	}else if(strcasecmp(identifier,"COLD") == 0) {
	  type = PPGPLOT_COLD;
	}else if(strcasecmp(identifier,"INVERTED_COLD") == 0) {
	  type = PPGPLOT_INVERTED_COLD;
	}else if(strcasecmp(identifier,"DIVREDBLUE") == 0) {
	  type = PPGPLOT_DIVREDBLUE;
	}else if(strcasecmp(identifier,"INVERTED_DIVREDBLUE") == 0 || strcasecmp(identifier,"NL") == 0) {
	  type = PPGPLOT_INVERTED_DIVREDBLUE;
//START REGION DEVELOP
	}else if(strcasecmp(identifier,"PLASMA") == 0) {
	  type = PPGPLOT_PLASMA;
	}else if(strcasecmp(identifier,"INVERTED_PLASMA") == 0) {
	  type = PPGPLOT_INVERTED_PLASMA;
	}else if(strcasecmp(identifier,"FOREST") == 0) {
	  type = PPGPLOT_FOREST;
	}else if(strcasecmp(identifier,"INVERTED_FOREST") == 0) {
	  type = PPGPLOT_INVERTED_FOREST;
	}else if(strcasecmp(identifier,"ALIEN_GLOW") == 0) {     // NOTE THAT THIS IS ALMOST THE SAME AS THE HEAT2 MAP
	  type = PPGPLOT_ALIEN_GLOW;
	}else if(strcasecmp(identifier,"INVERTED_ALIEN_GLOW") == 0) {
	  type = PPGPLOT_INVERTED_ALIEN_GLOW;
//START REGION RELEASE
	}else if(strcasecmp(identifier,"HEAT2") == 0) {
	  type = PPGPLOT_HEAT2;
	}else if(strcasecmp(identifier,"INVERTED_HEAT2") == 0) {
	  type = PPGPLOT_INVERTED_HEAT2;
	}else if(strcasecmp(identifier,"HEAT3") == 0) {
	  type = PPGPLOT_HEAT3;
	}else if(strcasecmp(identifier,"INVERTED_HEAT3") == 0) {
	  type = PPGPLOT_INVERTED_HEAT3;
	}else if(strcasecmp(identifier,"HEAT4") == 0) {
	  type = PPGPLOT_HEAT4;
	}else if(strcasecmp(identifier,"INVERTED_HEAT4") == 0) {
	  type = PPGPLOT_INVERTED_HEAT4;
	}else if(strcasecmp(identifier,"INFERNO") == 0 || strcasecmp(identifier,"CRISTINA") == 0) {
	  type = PPGPLOT_INFERNO;
	}else {
	  fflush(stdout);
	  printerror(debug, "ERROR cmap_parse_commandline:  '%s' not recognized as a color map.", identifier);
	  printCMAPCommandlineOptions(stderr);
	  return 0;	  
	}
      }
    }
  }
  return type;
}


//START REGION DEVELOP

/* 
   The angles dlat and dlong set the step size in latitude and
   longitude in degrees. lw sets the line width. Normally longitude 0
   and latitude 0 defines the centre of the projection, but this grid
   can be rotated by using rot_long and rot_lat, which are added to
   the longitude and latitude lines which are drawn.

   projection = 1 Draws a grid of a Hammer-Aitoff projection. 
   projection = 2 Draws a grid of a 3D sphere. 
   projection = 3 Draws a grid of a simple linear longitude-latitude map

   For projection 1 the poles are assumed to be at (0, +1) and (0, -1)
   and the equator extends from (-2, 0) to (+2, 0). For projection 2
   the sphere has a radius of 1. If a point lies on the front-side,
   the sphere is centred at x=-1.125, y=0. If the point lies on the
   far side the point ends on a sphere with an x-offset in the
   opposite direction. For projection 3 the horizontal axis runs from
   -180 to 180 degrees and the vertical axis from -90 to 90 degrees.
   */
void drawSphericalGrid(float dlat, float dlong, float rot_long, float rot_lat, int lw, int projection)
{
  float lat, lon, t, dt, x, y, sign, weight;
  int first, side, side2, oldside, oldside2, ok;

  ppgslw(lw);

  dlat *= M_PI/180.0;
  dlong *= M_PI/180.0;
  rot_long *= M_PI/180.0;
  rot_lat *= M_PI/180.0;

  /* The step size in radians */
  dt = 0.01;

  /* Draw outline */
  if(projection != 3) {
    for(side = 0; side < 2; side++) {
      first = 1;
      for(t = 0; t <= 2.0*M_PI+2.0*dt; t += dt) {
	ok = 1;
	if(projection == 1) {      /* Hammer-Aitoff projection: outline is a 2 by 1 ellipse */
	  x = 2*cos(t);
	  y = sin(t);
	  if(side == 1)
	    ok = 0;
	}else if(projection == 2) { /* Sphere projection: outline are two circles with radius 1 */
	  x = cos(t);
	  y = sin(t);
	  if(side == 0)
	    x -= 1.125;
	  else
	    x += 1.125;
	}
	if(first && ok) {
	  ppgmove(x, y);
	  first = 0;
	}else if(ok) {
	  ppgdraw(x, y);
	}
      }
    }
  }

  /* Draw latitude grid */
  for(sign = -1; sign < 1.1; sign += 2) {
    for(lat = 0; lat <= 0.5*M_PI+2.0*dt; lat += dlat) {
      first = 1;
      oldside = -1;
      oldside2 = -1;
      for(t = -M_PI; t <= M_PI+2.0*dt; t += dt) {
	lon = t;
	if(lon > M_PI)
	  lon = M_PI;
	if(projection == 1) {
	  projectionHammerAitoff_xy(lon, sign*lat, rot_long, rot_lat, &x, &y);
	  if(x < 0)
	    side = 0;
	  if(x > 0)
	    side = 1;
	  if(side != oldside)
	    first = 1;
	  oldside = side;
	}else if(projection == 2) {
	  side = projection_sphere_xy(lon, sign*lat, rot_long, rot_lat, &x, &y, &weight);
	  if(side != oldside)
	    first = 1;
	  oldside = side;
	}else if(projection == 3) {
	  projection_longlat_xy(lon, sign*lat, rot_long, rot_lat, &x, &y);
	  if(x < 0)
	    side = 0;
	  if(x > 0)
	    side = 1;
	  if(y < 0)
	    side2 = 0;
	  if(y > 0)
	    side2 = 1;
	  if(side != oldside || side2 != oldside2)
	    first = 1;
	  oldside = side;
	  oldside2 = side2;
	}
	if(first) {
	  ppgmove(x, y);
	  first = 0;
	}else {
	  ppgdraw(x, y);
	}
      }
    }
  }

  /* Draw longitude grid */
  for(sign = -1; sign < 1.1; sign += 2) {
    for(lon = 0; lon <= M_PI; lon += dlong) {
      first = 1;
      oldside = -1;
      oldside2 = -1;
      for(t = -0.5*M_PI; t <= 0.5*M_PI+2.0*dt; t += dt) {
	
	lat = t;
	if(lat > 0.5*M_PI)
	  lat = 0.5*M_PI;
	if(projection == 1) {
	  projectionHammerAitoff_xy(sign*lon, lat, rot_long, rot_lat, &x, &y);
	  if(x < 0)
	    side = 0;
	  if(x > 0)
	    side = 1;
	  if(side != oldside)
	    first = 1;
	  oldside = side;
	}else if(projection == 2) {
	  side = projection_sphere_xy(sign*lon, lat, rot_long, rot_lat, &x, &y, &weight);
	  if(side != oldside)
	    first = 1;
	  oldside = side;
	  //	  if(t >= 0) { printf("XXXXX side=%d lon=%f\n", side, sign*lon*180.0/M_PI); tmp = 0; }
	}else if(projection == 3) {
	  projection_longlat_xy(sign*lon, lat, rot_long, rot_lat, &x, &y);
	  if(x < 0)
	    side = 0;
	  if(x > 0)
	    side = 1;
	  if(y < 0)
	    side2 = 0;
	  if(y > 0)
	    side2 = 1;
	  if(side != oldside || side2 != oldside2)
	    first = 1;
	  oldside = side;
	  oldside2 = side2;
	}
	if(first) {
	  ppgmove(x, y);
	  first = 0;
	}else {
	  ppgdraw(x, y);
	}
      }
    }
  }

  ppgslw(1);
}

//START REGION DEVELOP
//START REGION RELEASE

/*
  Strips the device type from given device name, i.e. it returns:

  0  = undefined device type
  1  = /xwindow
  1  = /xw
  2  = /xserve
  2  = /xs
  3  = /ps
  4  = /vps
  5  = /cps
  6  = /vcps
  7  = /latex
  8  = /null
  9  = /png
  10 = /tpng
  100 = /pw

  If devicename = NULL, the name of the current opened device is used
  instead.

 */
int pgplot_device_type(char *devicename, verbose_definition verbose)
{
  int i, n, found, type;
  char *curdev;

  if(devicename != NULL) {  // If device name is supplied, use it
    curdev = devicename;
  }else {  // If device name not supplied, get it from pgplot
    curdev = calloc(1000, 1);
    if(curdev == NULL) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR pgplot_device_type: Memory allocation error");
      return 0;
    }
    n = 998;
    curdev[0] = '/';
    ppgqinf("type", curdev+1, &n);
    //    printf("device string length = %d\n", n);
    //    fflush(stdout);
    //    printf("Device = %s\n", curdev);
  }

  n = strlen(curdev);
  if(n <= 1)
    return 0;
  found = 0;
  for(i = n; i >= 0; i--) {
    if(curdev[i] == '/') {
      found = 1;
      break;
    }
  }


  type = 0;

  //printf("Device type = %s", &curdev[i]);

  if(found) {
    if(strcasecmp(&curdev[i], "/xwindow") == 0) {
      type = 1;
    }else if(strcasecmp(&curdev[i], "/xw") == 0) {
      type = 1;
    }else if(strcasecmp(&curdev[i], "/xserve") == 0) {
      type = 2; 
    }else if(strcasecmp(&curdev[i], "/xs") == 0) {
      type = 2;
    }else if(strcasecmp(&curdev[i], "/ps") == 0) {
      type = 3;
    }else if(strcasecmp(&curdev[i], "/vps") == 0) {
      type = 4;
    }else if(strcasecmp(&curdev[i], "/cps") == 0) {
      type = 5;
    }else if(strcasecmp(&curdev[i], "/vcps") == 0) {
      type = 6;
    }else if(strcasecmp(&curdev[i], "/latex") == 0) {
      type = 7;
    }else if(strcasecmp(&curdev[i], "/null") == 0) {
      type = 8;
    }else if(strcasecmp(&curdev[i], "/png") == 0) {
      type = 9;
    }else if(strcasecmp(&curdev[i], "/tpng") == 0) {
      type = 10;
    }else if(strcasecmp(&curdev[i], "/pw") == 0) {
      type = 100;
    }
  }

  //  printf(" (type %d)\n", type);

  if(devicename == NULL) {
    free(curdev);
  }
  return type;
}

//START REGION DEVELOP
