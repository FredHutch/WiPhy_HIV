//NOTE: This read-only file is generated from the file ../src/full_gui.cpp.
//Any editing should be done in that file.

/*** This file contains all the GUI code for real-time plotting using **/
/*** the stochastic early HIV model (with sequence and latent cell tracking)***/ 	 
#include<iostream>
#include<iomanip>
#include<fstream>
#include<string>
#include<cstdlib>
#include<cmath>
#include <map>
using namespace std;

#ifdef __sun
#include <strings.h>
#else
#include <string.h>
#endif

// Some STL Headers
#include <vector>
#include <stdlib.h>

// Using The STL Exception Library Increases The
// Chances That Someone Else Using Our Code Will Correctly
// Catch Any Exceptions That We Throw.
#include <stdexcept>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_blas.h>

#include <GL/gl.h>
#include <GL/glu.h>
#include <png.h>

#include <glib.h>
#include <gtk/gtk.h>
#include <gdk/gdkkeysyms.h>

#include <gtk/gtkgl.h>
#include <glib.h>
#include <pango/pangocairo.h>

#include "colors.h"  
#include "plotpoints.h"  
#include "settings.h"  
#include "strain.h"  
#include "gui.h"  

#include "hiv.xpm"

static GLfloat Black[4] = {0,0,0,1};
static GLfloat Pink[4] = {1, 0.5, 0.5, 1};
static GLfloat DarkBlue[4] = {0., 0., 0.65, 1};
static GLfloat DarkRed[4] = {0.65, 0., 0., 1};
static GLfloat DarkGreen[4] = {0.0, 0.65, 0.0, 1};
static GLfloat DarkPurple[4] = {0.44, 0., 0.73, 1};
static GLfloat Yellow[4] = {1, 1, 0.2, 1};
static GLfloat White[4] = {1,1,1,1};
static GLfloat Gray[4] = {.75,.75,1,1};
static GLfloat Orange[4] = {0.75, 0.375, 0.0, 1};
static GLfloat Blue[4] = {0.25, 0.25, 1.00, 1};
static GLfloat Gold[4] = {1.0, 211./255., 32./255., 1};
static GLfloat Brick[4] = {126./255., 0.0, 33./255., 1};

using namespace std;

#define MIN(a,b) (((a) < (b)) ? (a) : (b))
#define DEFAULT_WIDTH  800
#define DEFAULT_HEIGHT 600
#define DEFAULT_IMAGE_WIDTH  816
#define DEFAULT_IMAGE_HEIGHT 592
#ifdef __sun
#define DEFAULT_TITLE  "HIV simulation (for SunOS)"
#else
#define DEFAULT_TITLE  "HIV simulation"
#endif
#define PROGRAM_NAME  "hiv_sim"
#define PACKAGE_VERSION  "0.1"
#define PACKAGE_BUGREPORT  "dswan@fredhutch.org"

#define ORIGIN_X          0
#define ORIGIN_Y          0
#define MAX_X_COORD          80.
#define MAX_Y_COORD          80.

#define MAX_RGB_VAL 255
#define MAX_GUI_ENTRIES 50

#define VT_GRAPH 1
#define CD8_GRAPH 2
#define ACT_GRAPH 3
#define SHOW_PARAMS 4
#define SHOW_STATS 5

static settings *theState=NULL;
static int drawnTime = 0;

static void *bgr;
static bool draw_routine (GLfloat width, GLfloat height);
void write_png_file( int width, int height,const char* file_name);
static string outDir;

static GLuint font_list_base;

static gchar default_font_string[] = "Courier Bold 14";

static GLuint default_font_list_base;

static PangoFontDescription *default_font_desc = NULL;

static GdkGLContext *global_glcontext;
static GdkGLDrawable *global_gldrawable;

static GtkWidget *main_window;

static GtkWidget *image;
static GdkGLConfig *glconfig;

static GtkWidget *ART_start_entry;
static GtkWidget *ART_stop_entry;
static GtkWidget *ART_eff_entry;
static GtkWidget *beta_entry;
static GtkWidget *beta_max_entry;
static GtkWidget *beta_k_entry;
static GtkWidget *lifespan_decay_entry;
static GtkWidget *gam_entry;
static GtkWidget *mu_entry;
static GtkWidget *pi_entry;
static GtkWidget *theta_entry;
static GtkWidget *junk_entry;
static GtkWidget *dA_ic50_entry;
static GtkWidget *dA_scale_entry;
static GtkWidget *eclipse_entry;
static GtkWidget *delta_entry;
static GtkWidget *cd8_mean_entry;
static GtkWidget *cd8_0_entry;
static GtkWidget *cd8_k_entry;
static GtkWidget *dImmun_entry;
static GtkWidget *dImmun_s_entry;
static GtkWidget *dImmun_IC50_entry;
static GtkWidget *kappa_entry;
static GtkWidget *cd8_ic50_entry;
static GtkWidget *Immun_model_entry;
static GtkWidget *followTopStrains;
static GtkWidget *plotTopStrains;
static GtkWidget *plotFirstStrains;
static GtkWidget *colorByHamming;
static GtkWidget *maxHamming_entry;
static GtkWidget *deltaHamming_entry;

static GtkWidget *backup_slider;
static GtkObject *adjuster;

static GtkWidget *font_name_widget;
static GtkWidget *font_size_widget;

static GtkWidget *background_r;
static GtkWidget *background_g;
static GtkWidget *background_b;
static GtkWidget *rgb_box;

static GtkWidget *zoomin_button;
static GtkWidget *zoomout_button;
static GtkWidget *run_button;
static GtkWidget *stop_button;
static GtkWidget *pause_button;
static GtkWidget *page_left_button;
static GtkWidget *page_right_button;
static GtkWidget *strain_widget;
static GtkWidget *cell_widget;
static GtkWidget *top_strain_box;
static GtkWidget *first_strain_box;

static GtkWidget *max_log_v_entry;
static GtkWidget *max_log_a_entry;
static GtkWidget *max_log_cd8_entry;

static GLfloat da_width = 0.0;
static GLfloat da_height = 0.0;

static	int enc_width, enc_height;
static bool onePrint = false;
static bool legend = false;


/**************************************************************************
 *  * The following section contains the function prototype declarations.
 *   **************************************************************************/

void updateGL(void);

static void stop_cb (GtkWidget* , gpointer );

void wait_cursor( GdkWindow *win)
{
    GdkCursor *cur;
    cur = gdk_cursor_new( GDK_CLOCK );
    gdk_window_set_cursor( win, cur );
    gdk_cursor_unref( cur );
}

void normal_cursor( GdkWindow *win)
{
    gdk_window_set_cursor( win, NULL );
}

void abort_(const char * s, ...)
{
        va_list args;
        va_start(args, s);
        vfprintf(stderr, s, args);
        fprintf(stderr, "\n");
        va_end(args);
        abort();
}
static GLfloat* pColors[MAX_STRAIN_COLORS+2];
static GLfloat* rainbow[MAX_RAINBOW_COLORS];

static int simLock = 0;

void assignRainbow(GLfloat *color_arr, int index, int max_index)
{
    GLfloat red;
    GLfloat green;
    GLfloat blue;

    if (index < (double)max_index/3) {
	red = 1;
	green=((index/((double)max_index/3))*1.0);
	blue=0;
    } else if (index < max_index/2) {
	red=(((double)max_index/2 - index)/((double)max_index/6))*1.0;
	green = 1.0;
	blue=0;
    } else if (index < 2*max_index/3) {
	red=0;
	blue=((index-(double)max_index/2)/((double)max_index/6))*1.0;
	green = 1.0;
    } else {
	red = 0;
	green=(((double)(max_index - index))/((double)max_index/3))*1.0;
	blue=1.0;
    }
    color_arr[0] = red;
    color_arr[1] = green;
    color_arr[2] = blue;
    color_arr[3] = 1.0;
}

void take_screenshot( GtkWidget *ok_button, char *filename );
double ScoreFunction(settings *vars);

/* This function runs one simulation.  It is the entry point when 
 * the program is run in the non-batch mode */
void runSimulation ( settings *vars)
{
    ScoreFunction(vars);
}

/**************************************************************************
 * The following section contains all the callback function definitions.
 *  **************************************************************************/

png_byte color_type;
png_byte bit_depth;

png_structp png_ptr;
png_infop info_ptr;
int number_of_passes;
png_bytep * row_pointers;

void write_png_file( int width, int height,const char* file_name)
{
        /* create file */
        FILE *fp = fopen(file_name, "wb");
        if (!fp)
                abort_("[write_png_file] File %s could not be opened for writing", file_name);

        /* initialize stuff */
        png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);

        if (!png_ptr)
                abort_("[write_png_file] png_create_write_struct failed");

        info_ptr = png_create_info_struct(png_ptr);
        if (!info_ptr)
                abort_("[write_png_file] png_create_info_struct failed");

	row_pointers = (png_bytep*) malloc(sizeof(png_bytep) * height);

	char *currentP=(char *)bgr;
	for (int j=0; j < height; j++)
	{
	    row_pointers[height-j-1]=(png_bytep)currentP;
	    currentP+=width*4;
	}

        if (setjmp(png_jmpbuf(png_ptr)))
                abort_("[write_png_file] Error during init_io");

        png_init_io(png_ptr, fp);


        /* write header */
        if (setjmp(png_jmpbuf(png_ptr)))
                abort_("[write_png_file] Error during writing header");

        png_set_IHDR(png_ptr, info_ptr, width, height,
                     8, PNG_COLOR_TYPE_RGB_ALPHA, PNG_INTERLACE_NONE,
                     PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);

        png_write_info(png_ptr, info_ptr);


        /* write bytes */
        if (setjmp(png_jmpbuf(png_ptr)))
                abort_("[write_png_file] Error during writing bytes");

        png_write_image(png_ptr, row_pointers);


        /* end write */
        if (setjmp(png_jmpbuf(png_ptr)))
                abort_("[write_png_file] Error during end of write");

        png_write_end(png_ptr, NULL);

        /* cleanup heap allocation */
        free(row_pointers);

        fclose(fp);
}
#ifndef NO_GUI

bool ShowMessageBox(std::string , std::string message)
{
    GtkWidget*dialog = gtk_message_dialog_new (GTK_WINDOW(main_window),
                                 GTK_DIALOG_DESTROY_WITH_PARENT,
                                 GTK_MESSAGE_ERROR,
                                 GTK_BUTTONS_CLOSE,
                                 message.c_str());

    gtk_dialog_run (GTK_DIALOG (dialog));
    gtk_widget_destroy (dialog);

    return true;
}

void take_screenshot( GtkWidget *, char *filename )
{
	GdkPixbuf *screenshot = NULL;

	screenshot = gdk_pixbuf_get_from_drawable( screenshot, GDK_DRAWABLE(GTK_WIDGET(image)->window), gdk_colormap_get_system(), 0, 0, 0, 0, enc_width, enc_height );
	gdk_pixbuf_save( GDK_PIXBUF(screenshot), filename, "png", NULL, NULL );
	//g_free( filename );
}


void file_name_changed( GtkWidget *widget, GtkWidget *check_button )
{
	if( strlen( gtk_entry_get_text(GTK_ENTRY(widget)) ) < 1 )
		gtk_widget_set_sensitive( GTK_WIDGET(check_button), FALSE );
	else
		gtk_widget_set_sensitive( GTK_WIDGET(check_button), TRUE );
}
#ifdef SHOW_FONT_LIST
static void
list_fonts ()
{
    int i;
    PangoFontFamily ** families;
    int n_families;
    PangoFontMap * fontmap;

    fontmap = pango_cairo_font_map_get_default();
    pango_font_map_list_families (fontmap, & families, & n_families);
    printf ("There are %d families\n", n_families);
    for (i = 0; i < n_families; i++) {
        PangoFontFamily * family = families[i];
        const char * family_name;

        family_name = pango_font_family_get_name (family);
        printf ("Family %d: %s\n", i, family_name);
    }
    g_free (families);
}
#endif

void font_changed( GtkWidget *)
{
    PangoFontDescription *font_desc;
    PangoFont *font;
    string font_name;
    string font_size;


    font_name = gtk_entry_get_text(GTK_ENTRY(font_name_widget));
    font_size = gtk_entry_get_text(GTK_ENTRY(font_size_widget));

    font_name += " "+font_size;

    if( font_name.length() > 0 )
    {
	font_desc = pango_font_description_from_string (font_name.c_str());
        font_list_base = glGenLists (128);

	font = gdk_gl_font_use_pango_font (font_desc, 0, 128, font_list_base);
	if (font == NULL)
	{
	    string err_str ="*** Can't load font '"+font_name+"'\n";
	    ShowMessageBox("Warning",err_str);

	    g_print (err_str.c_str());
	    glDeleteLists(font_list_base,128);

	    font_list_base = default_font_list_base;
	    font_desc = default_font_desc;
	    font = gdk_gl_font_use_pango_font (font_desc, 0, 128, font_list_base);
	}
	else
	{
	    g_print ("*** Loaded font '%s'\n", font_name.c_str());
	    glDeleteLists(default_font_list_base,128);
	    default_font_list_base = font_list_base;
	    default_font_desc = font_desc;
	}
    }
}
void background_changed(GtkWidget *)
{
    double back_r = (double)gtk_spin_button_get_value_as_int( GTK_SPIN_BUTTON(background_r) );
    double back_g = (double)gtk_spin_button_get_value_as_int( GTK_SPIN_BUTTON(background_g) );
    double back_b = (double)gtk_spin_button_get_value_as_int( GTK_SPIN_BUTTON(background_b) );

    glClearColor(back_r/255., back_g/255., back_b/255., 1.0);

    GdkColor init_color;

    init_color.red = 65535 * (back_r/255.);;
    init_color.green = 65535 * (back_g/255.);
    init_color.blue = 65535 * (back_b/255.);

    gtk_widget_modify_bg(rgb_box,GTK_STATE_NORMAL, &init_color);
}


void open_file_selection( GtkWidget *widget, gpointer data )
{
	GtkWidget *dialog;
	char *filename = NULL;
	void (*filename_handler)(GtkWidget *, char *) = (void (*)(GtkWidget*, char*))data;
	/* Function to call if a filename selected */

	dialog = gtk_file_chooser_dialog_new("Select File...", NULL, GTK_FILE_CHOOSER_ACTION_SAVE, GTK_STOCK_CANCEL, GTK_RESPONSE_CANCEL, GTK_STOCK_OPEN, GTK_RESPONSE_ACCEPT, NULL);
	if( gtk_dialog_run(GTK_DIALOG(dialog)) == GTK_RESPONSE_ACCEPT ) {
		filename = gtk_file_chooser_get_filename( GTK_FILE_CHOOSER(dialog) );
	}
	gtk_widget_destroy( dialog );
	if( filename != NULL )
		filename_handler( widget, filename );
}


static gboolean
initialize_image(GLfloat width, GLfloat height)
{
	//gtk_widget_set_size_request( GTK_WIDGET(image), width, height );
	global_gldrawable = gtk_widget_get_gl_drawable( image );
	global_glcontext = gtk_widget_get_gl_context( image );

	if( gdk_gl_drawable_gl_begin(global_gldrawable, global_glcontext) ) {

		glClear( GL_COLOR_BUFFER_BIT );
		glViewport(0, 0, width, height);
		glMatrixMode( GL_PROJECTION );
		 glDisable ( GL_LIGHTING ) ;
		glLoadIdentity();
		glOrtho(-60.0, 60.0,-60., 60.,  1.0, -1.0);
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
		glFlush();

		if (gdk_gl_drawable_is_double_buffered (global_gldrawable))
		  gdk_gl_drawable_swap_buffers (global_gldrawable);
		gdk_gl_drawable_gl_end( global_gldrawable );
	}
	gtk_widget_queue_draw( image );

	/* Encoders encode better with heights or widths that are multiples of 16 */
	/* so the area to be encoded is clipped to satisfy this. */
	enc_width = (int)width - (int)width%16;
	enc_height = (int)height - (int)height%16;
	return true;
}
void select_font(GtkWidget *widget, gpointer label);
/***
 *  *** The "realize" signal handler. All the OpenGL initialization
 *   *** should be performed here, such as default background colour,
 *    *** certain states etc.
 *     ***/
static void
realize (GtkWidget *widget,
         gpointer   )
{
  GdkGLContext *glcontext = gtk_widget_get_gl_context (widget);
  GdkGLDrawable *gldrawable = gtk_widget_get_gl_drawable (widget);

  GLfloat width = widget->allocation.width;
  GLfloat height = widget->allocation.height;

  da_width = widget->allocation.width;
  da_height = widget->allocation.height;

  PangoFontDescription *font_desc;
  PangoFont *font;

  /*** OpenGL BEGIN ***/
  if (!gdk_gl_drawable_gl_begin (gldrawable, glcontext))
    return;

  /*
   * Generate font display lists.
   */
  if (default_font_desc == NULL)
  {
        default_font_desc = pango_font_description_from_string (default_font_string);
        default_font_list_base = glGenLists (128);
  }

  font_list_base = default_font_list_base;

  font_desc = default_font_desc;

  font = gdk_gl_font_use_pango_font (font_desc, 0, 128, font_list_base);
  if (font == NULL)
    {
      g_print ("*** Can't load font '%s'\n", default_font_string);
      select_font(NULL, NULL);
      //exit (1);
    }

  /*** Fill in the details here. ***/
  //glClearColor(170./255., 180./255., 190./255., 1.0);

  double back_r = (double)gtk_spin_button_get_value_as_int( GTK_SPIN_BUTTON(background_r) );
  double back_g = (double)gtk_spin_button_get_value_as_int( GTK_SPIN_BUTTON(background_g) );
  double back_b = (double)gtk_spin_button_get_value_as_int( GTK_SPIN_BUTTON(background_b) );

  glClearColor(back_r/255., back_g/255., back_b/255., 1.0);
  glClearDepth(1.0);
  /*glShadeModel(GL_FLAT);*/
  glDisable(GL_DEPTH_TEST);

  gdk_gl_drawable_gl_end (gldrawable);
  /*** OpenGL END ***/

  initialize_image(width, height);

  return;
}

static gboolean
configure_event (GtkWidget         *widget,
                 GdkEventConfigure *)
{
	GLfloat width = widget->allocation.width;
	GLfloat height = widget->allocation.height;
	gboolean value = TRUE;

	if (width != da_width || height != da_height)
		realize(widget,NULL);
	else
		value =  initialize_image(width, height);

	return value;
}
void setCheckVt(GtkWidget *widget, gpointer )
{
    int value = (gtk_toggle_button_get_active (GTK_TOGGLE_BUTTON(widget)))?1:0;
    if (theState != NULL) {
	theState->plotVt = value;
        updateGL();
    }
}

void setCheckColorByHamming(GtkWidget *widget, gpointer )
{
    int value = (gtk_toggle_button_get_active (GTK_TOGGLE_BUTTON(widget)))?1:0;
    if (theState != NULL) {
	theState->colorByHamming = value;
        updateGL();
    }
}

void setCheckAct(GtkWidget *widget, gpointer )
{
    int value = (gtk_toggle_button_get_active (GTK_TOGGLE_BUTTON(widget)))?1:0;
    if (theState != NULL) {
	theState->plotAct = value;
        updateGL();
    }
}

void setCheckCD8s(GtkWidget *widget, gpointer )
{
    int value = (gtk_toggle_button_get_active (GTK_TOGGLE_BUTTON(widget)))?1:0;
    if (theState != NULL) {
	theState->plotCD8s = value;
        updateGL();
    }
}

void setShowParams(GtkWidget *widget, gpointer )
{
    int value = (gtk_toggle_button_get_active (GTK_TOGGLE_BUTTON(widget)))?1:0;
    if (theState != NULL) {
	theState->show_params = value;
        updateGL();
    }
}

void setShowStats(GtkWidget *widget, gpointer )
{
    int value = (gtk_toggle_button_get_active (GTK_TOGGLE_BUTTON(widget)))?1:0;
    if (theState != NULL) {
	theState->show_stats = value;
        updateGL();
    }
}

static void
setMaxStrains (GtkWidget* widget, gpointer )
{
    int value;
    char value_str[100];
    strcpy(value_str,gtk_entry_get_text(GTK_ENTRY(widget)));
    sscanf(value_str,"%d", &value);

    if (theState != NULL) {
	theState->maxTopStrains = value;
    }
    updateGL();
}
static void
updateParams (GtkWidget* widget, gpointer )
{
    char value_str[100];
    strcpy(value_str,gtk_entry_get_text (GTK_ENTRY(ART_eff_entry)));
    sscanf(value_str,"%lf", &theState->ART_eff);
    strcpy(value_str,gtk_entry_get_text (GTK_ENTRY(ART_start_entry)));
    sscanf(value_str,"%lf", &theState->ART_start);
    strcpy(value_str,gtk_entry_get_text (GTK_ENTRY(ART_stop_entry)));
    sscanf(value_str,"%lf", &theState->ART_stop);
    strcpy(value_str,gtk_entry_get_text (GTK_ENTRY(beta_entry)));
    sscanf(value_str,"%lf", &theState->Bt_mean);
    strcpy(value_str,gtk_entry_get_text (GTK_ENTRY(beta_max_entry)));
    sscanf(value_str,"%lf", &theState->Bt_max);
    strcpy(value_str,gtk_entry_get_text (GTK_ENTRY(beta_k_entry)));
    sscanf(value_str,"%lf", &theState->Bt_k);
    strcpy(value_str,gtk_entry_get_text (GTK_ENTRY(lifespan_decay_entry)));
    sscanf(value_str,"%lf", &theState->lifespan_decay);
    strcpy(value_str,gtk_entry_get_text (GTK_ENTRY(gam_entry)));
    sscanf(value_str,"%lf", &theState->gam);
    strcpy(value_str,gtk_entry_get_text (GTK_ENTRY(junk_entry)));
    sscanf(value_str,"%lf", &theState->junk);
    strcpy(value_str,gtk_entry_get_text (GTK_ENTRY(mu_entry)));
    sscanf(value_str,"%lf", &theState->mu0);
    strcpy(value_str,gtk_entry_get_text (GTK_ENTRY(pi_entry)));
    sscanf(value_str,"%lf", &theState->pi);
    strcpy(value_str,gtk_entry_get_text (GTK_ENTRY(theta_entry)));
    sscanf(value_str,"%lf", &theState->theta);
    strcpy(value_str,gtk_entry_get_text (GTK_ENTRY(dA_ic50_entry)));
    sscanf(value_str,"%lf", &theState->dA_ic50_0);
    strcpy(value_str,gtk_entry_get_text (GTK_ENTRY(dA_scale_entry)));
    sscanf(value_str,"%lf", &theState->dA_scale);
    strcpy(value_str,gtk_entry_get_text (GTK_ENTRY(eclipse_entry)));
    sscanf(value_str,"%lf", &theState->eclipse);
    strcpy(value_str,gtk_entry_get_text (GTK_ENTRY(delta_entry)));
    sscanf(value_str,"%lf", &theState->delta);
    strcpy(value_str,gtk_entry_get_text (GTK_ENTRY(cd8_mean_entry)));
    sscanf(value_str,"%lf", &theState->cd8_mean);
    strcpy(value_str,gtk_entry_get_text (GTK_ENTRY(cd8_k_entry)));
    sscanf(value_str,"%lf", &theState->cd8_k);
    strcpy(value_str,gtk_entry_get_text (GTK_ENTRY(cd8_0_entry)));
    sscanf(value_str,"%lf", &theState->cd8_0);
    strcpy(value_str,gtk_entry_get_text (GTK_ENTRY(dImmun_entry)));
    sscanf(value_str,"%lf", &theState->dImmun);
    strcpy(value_str,gtk_entry_get_text (GTK_ENTRY(dImmun_IC50_entry)));
    sscanf(value_str,"%lf", &theState->dImmun_IC50_0);
    strcpy(value_str,gtk_entry_get_text (GTK_ENTRY(dImmun_s_entry)));
    sscanf(value_str,"%lf", &theState->dImmun_s);
    strcpy(value_str,gtk_entry_get_text (GTK_ENTRY(kappa_entry)));
    sscanf(value_str,"%lf", &theState->kappa);
    strcpy(value_str,gtk_entry_get_text (GTK_ENTRY(cd8_ic50_entry)));
    sscanf(value_str,"%lf", &theState->cd8_ic50_0);
    strcpy(value_str,gtk_entry_get_text (GTK_ENTRY(dA_ic50_entry)));
    sscanf(value_str,"%lf", &theState->dA_ic50_0);
    strcpy(value_str,gtk_entry_get_text (GTK_ENTRY(eclipse_entry)));
    sscanf(value_str,"%lf", &theState->eclipse);
    strcpy(value_str,gtk_entry_get_text (GTK_ENTRY(Immun_model_entry)));
    sscanf(value_str,"%d", &theState->immun_model);
    strcpy(value_str,gtk_entry_get_text (GTK_ENTRY(strain_widget)));
    sscanf(value_str,"%d", &theState->maxTopStrains);
    strcpy(value_str,gtk_entry_get_text (GTK_ENTRY(cell_widget)));
    sscanf(value_str,"%d", &theState->cellThreshold);
    strcpy(value_str,gtk_entry_get_text (GTK_ENTRY(max_log_v_entry)));
    sscanf(value_str,"%lf", &theState->max_log_v_value);
    strcpy(value_str,gtk_entry_get_text (GTK_ENTRY(max_log_a_entry)));
    sscanf(value_str,"%lf", &theState->max_log_a_value);
    strcpy(value_str,gtk_entry_get_text (GTK_ENTRY(max_log_cd8_entry)));
    sscanf(value_str,"%lf", &theState->max_log_cd8_value);
    strcpy(value_str,gtk_entry_get_text (GTK_ENTRY(maxHamming_entry)));
    sscanf(value_str,"%d", &theState->maxHammingDist);
    strcpy(value_str,gtk_entry_get_text (GTK_ENTRY(deltaHamming_entry)));
    sscanf(value_str,"%d", &theState->cd8GroupSize);
    theState->colorByHamming = (gtk_toggle_button_get_active (GTK_TOGGLE_BUTTON(colorByHamming)))?1:0;

    updateGL();
}

static void
setCellThreshold (GtkWidget* widget, gpointer )
{
    int value;
    char value_str[100];
    strcpy(value_str,gtk_entry_get_text(GTK_ENTRY(widget)));
    sscanf(value_str,"%d", &value);

    if (theState != NULL) {
	theState->cellThreshold = value;
    }
    updateGL();
}

void check_button_clicked(GtkWidget *widget, gpointer )
{
    int value = (gtk_toggle_button_get_active (GTK_TOGGLE_BUTTON(widget)))?1:0;
    if (theState != NULL) {
	theState->AutoSnapshot = value;
        updateGL();
    }
}

void setCheckFollowStrains(GtkWidget *widget, gpointer )
{
    int value = (gtk_toggle_button_get_active (GTK_TOGGLE_BUTTON(widget)))?1:0;
    if (theState != NULL) {
	theState->followTopStrains = value;
	if (theState->followTopStrains == 0 && theState->plotFirstStrains == 0 && theState->plotTopStrains == 0)
	{
	    gtk_widget_hide( GTK_WIDGET(top_strain_box) );
	    gtk_widget_hide( GTK_WIDGET(first_strain_box) );
	}
	else
	{
	    gtk_widget_show( GTK_WIDGET(top_strain_box) );
	    if (theState->plotFirstStrains == 0)
		gtk_widget_hide( GTK_WIDGET(first_strain_box) );
	    else
		gtk_widget_show( GTK_WIDGET(first_strain_box) );
	}
	if (theState->followTopStrains)
	{
	    theState->plotFirstStrains = 0;
	    theState->plotTopStrains = 0;
	}

	gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (plotFirstStrains), theState->plotFirstStrains != 0);
	gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (plotTopStrains), theState->plotTopStrains != 0);
        updateGL();
    }
}

void setCheckTopStrains(GtkWidget *widget, gpointer )
{
    int value = (gtk_toggle_button_get_active (GTK_TOGGLE_BUTTON(widget)))?1:0;
    if (theState != NULL) {
	theState->plotTopStrains = value;
	if (theState->followTopStrains == 0 && theState->plotFirstStrains == 0 && theState->plotTopStrains == 0)
	{
	    gtk_widget_hide( GTK_WIDGET(top_strain_box) );
	}
	else
	{
	    gtk_widget_show( GTK_WIDGET(top_strain_box) );
	}
	if (theState->plotTopStrains)
	{
	    theState->plotFirstStrains = 0;
	    theState->followTopStrains = 0;
	}

	if (theState->plotFirstStrains == 0)
	    gtk_widget_hide( GTK_WIDGET(first_strain_box) );
	else
	    gtk_widget_show( GTK_WIDGET(first_strain_box) );

	gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (plotFirstStrains), theState->plotFirstStrains != 0);
	gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (followTopStrains), theState->followTopStrains != 0);
        updateGL();
    }
}
void setCheckFirstStrains(GtkWidget *widget, gpointer )
{
    int value = (gtk_toggle_button_get_active (GTK_TOGGLE_BUTTON(widget)))?1:0;
    if (theState != NULL) {
	theState->plotFirstStrains = value;
	if (theState->followTopStrains == 0 && theState->plotFirstStrains == 0 && theState->plotTopStrains == 0)
	{
	    gtk_widget_hide( GTK_WIDGET(first_strain_box) );
	    gtk_widget_hide( GTK_WIDGET(top_strain_box) );
	}
	else
	{
	    gtk_widget_show( GTK_WIDGET(top_strain_box) );
	}
	if (theState->plotFirstStrains)
	{
	    theState->plotTopStrains = 0;
	    theState->followTopStrains = 0;
	    gtk_widget_show( GTK_WIDGET(first_strain_box) );
	}
	else
	    gtk_widget_hide( GTK_WIDGET(first_strain_box) );

	gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (plotTopStrains), theState->plotTopStrains != 0);
	gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (followTopStrains), theState->followTopStrains != 0);
        updateGL();
    }
}

void setCheckWriteOn(GtkWidget *widget, gpointer )
{
    int value = (gtk_toggle_button_get_active (GTK_TOGGLE_BUTTON(widget)))?1:0;
    if (theState != NULL) {
	theState->writeOn = value;
        updateGL();
    }
    onePrint=true;
}

void setCheckScrollAxes(GtkWidget *widget, gpointer )
{
    int value = (gtk_toggle_button_get_active (GTK_TOGGLE_BUTTON(widget)))?1:0;
    if (theState != NULL) {
	theState->scrollAxes = value;
        updateGL();
    }
}

static double plot_bias_saved=0.;
static void
run_cb (GtkWidget* , gpointer )
{
    gtk_widget_set_sensitive( GTK_WIDGET(run_button), FALSE );
    gtk_widget_set_sensitive( GTK_WIDGET(pause_button), TRUE );
    gtk_widget_set_sensitive( GTK_WIDGET(stop_button), TRUE );

    if (theState->pauseFlag == 1) {
	theState->pauseFlag = 0;
	theState->plot_bias = plot_bias_saved;
	gtk_widget_set_sensitive( GTK_WIDGET(page_left_button), FALSE );
	gtk_widget_set_sensitive( GTK_WIDGET(page_right_button), FALSE );
	gtk_widget_set_sensitive( GTK_WIDGET(followTopStrains), FALSE );
	gtk_widget_set_sensitive( GTK_WIDGET(plotTopStrains), FALSE );
	gtk_widget_hide( GTK_WIDGET(backup_slider) );
	updateGL();
	return;
    }

    // try to prevent double clicking or premature sim launches
    if (simLock != 0) {
	return;
    }


    simLock = 1;

    if (theState->points != NULL)
	delete theState->points;

    theState->points = new plotPoints();

    theState->stopFlag = 0;

    theState->plot_bias = 0.;
    theState->time_bias = 0;

    gtk_widget_set_sensitive( GTK_WIDGET(followTopStrains), FALSE );
    gtk_widget_set_sensitive( GTK_WIDGET(plotTopStrains), FALSE );
    gtk_widget_set_sensitive( GTK_WIDGET(plotFirstStrains), FALSE );
    gtk_widget_set_sensitive( GTK_WIDGET(strain_widget), FALSE );
    gtk_widget_set_sensitive( GTK_WIDGET(cell_widget), FALSE );

    /* reset scroll bar if necessary */
    ((GtkAdjustment *)adjuster)->value = 20;
    gtk_signal_emit_by_name(GTK_OBJECT(adjuster), "value_changed");

    gtk_widget_hide( GTK_WIDGET(backup_slider) );

    updateParams(NULL,NULL);

    runSimulation(theState);  

    if(theState->stopFlag != 0)
	stop_cb (NULL, NULL);
    gtk_widget_set_sensitive( GTK_WIDGET(followTopStrains), TRUE );
    gtk_widget_set_sensitive( GTK_WIDGET(plotTopStrains), TRUE );
    gtk_widget_set_sensitive( GTK_WIDGET(plotFirstStrains), TRUE );
    gtk_widget_set_sensitive( GTK_WIDGET(strain_widget), TRUE );
    gtk_widget_set_sensitive( GTK_WIDGET(cell_widget), TRUE );

    simLock = 0;
    updateGL();
}
static void
pause_cb (GtkWidget* , gpointer )
{
    // disallow double clicking or premature sim stopping
    if (theState->pauseFlag ==1) {
	//theState->app->beep();
	return;
    }

    theState->pauseFlag = 1;
    plot_bias_saved = theState->plot_bias;
    //gtk_widget_set_sensitive( GTK_WIDGET(zoomin_button), TRUE );
    //gtk_widget_set_sensitive( GTK_WIDGET(zoomout_button), TRUE );
    gtk_widget_set_sensitive( GTK_WIDGET(page_left_button), TRUE );
    gtk_widget_set_sensitive( GTK_WIDGET(page_right_button), TRUE );
    gtk_widget_set_sensitive( GTK_WIDGET(pause_button), FALSE );
    gtk_widget_set_sensitive( GTK_WIDGET(run_button), TRUE );
    gtk_widget_show( GTK_WIDGET(backup_slider) );
    updateGL();
}

static void
stop_cb (GtkWidget* , gpointer )
{
    // disallow double clicking or premature sim stopping
    if (simLock == 0) {
	//theState->app->beep();
	return;
    }

    theState->stopFlag = 1;
    theState->pauseFlag = 0;

    //gtk_widget_set_sensitive( GTK_WIDGET(zoomin_button), TRUE );
    //gtk_widget_set_sensitive( GTK_WIDGET(zoomout_button), TRUE );
    gtk_widget_set_sensitive( GTK_WIDGET(page_left_button), TRUE );
    gtk_widget_set_sensitive( GTK_WIDGET(page_right_button), TRUE );
    gtk_widget_set_sensitive( GTK_WIDGET(stop_button), FALSE );
    gtk_widget_set_sensitive( GTK_WIDGET(pause_button), FALSE );
    gtk_widget_set_sensitive( GTK_WIDGET(run_button), TRUE );
    gtk_widget_set_sensitive( GTK_WIDGET(followTopStrains), TRUE );
    gtk_widget_set_sensitive( GTK_WIDGET(plotTopStrains), TRUE );
    gtk_widget_set_sensitive( GTK_WIDGET(plotFirstStrains), TRUE );
    gtk_widget_hide( GTK_WIDGET(backup_slider) );
    updateGL();
}

static void
zoomin_cb (GtkWidget* , gpointer )
{
    if (theState->pauseFlag ==0) {
	//theState->app->beep();
	//return;
    }
    theState->plot_span = theState->plot_span/2.0;

    //emit param8Changed(QString::number(theState->plot_span));
    updateGL();
}
static void
zoomout_cb (GtkWidget* , gpointer )
{
    if (theState->pauseFlag ==0) {
	//theState->app->beep();
	//return;
    }
    theState->plot_span = theState->plot_span*2.0;

    //emit param8Changed(QString::number(theState->plot_span));
    updateGL();
}

static void
page_left_cb (GtkWidget* , gpointer )
{
    // disallow scrolling during run mode
    if (theState->pauseFlag !=1 && theState->stopFlag != 1) {
	//theState->app->beep();
	return;
    }
    if (theState->time+(theState->plot_bias- theState->plot_span) >= 0.0)
	theState->plot_bias -= theState->plot_span;
    updateGL();
}
static void
page_right_cb (GtkWidget* , gpointer )
{
    // disallow scrolling during run mode
    if (theState->pauseFlag !=1 && theState->stopFlag != 1) {
	//theState->app->beep();
	return;
    }
    if (theState->plot_bias < 0.0)
	theState->plot_bias += theState->plot_span;

    updateGL();
}

void scrollTime(GtkWidget* , gpointer )
{
    double value;

    // disallow scrolling during run mode
    if (theState->pauseFlag !=1 && theState->stopFlag != 1) {
	//theState->app->beep();
	return;
    }
    value = gtk_adjustment_get_value(GTK_ADJUSTMENT(adjuster));

    theState->time_bias = value - 20.0;

    if (theState->time_bias < 0 && theState->time+theState->time_bias < 0)
	theState->time_bias = -theState->time;

    updateGL();
}

void updateGL(void)
{
    //advance_handler();
    gdk_window_invalidate_rect (image->window, &image->allocation, FALSE);
}
//! [8-]
void resizeGL(int width, int height)
{
//    int side = qMin(width, height);
//    glViewport((width - side) / 2, (height - side) / 2, side, side);

    float aspectRatio=(float)width/(float)height;

    glViewport(0, 0, width, height);
    glMatrixMode(GL_PROJECTION);
		 glDisable ( GL_LIGHTING ) ;
    glLoadIdentity();
if (width <= height)
    glOrtho(-60., 60., -60.0/aspectRatio, 60.0/aspectRatio, 1.0, -1.0);
else
    glOrtho(-60.0*aspectRatio, 60.0*aspectRatio,-60., 60.,  1.0, -1.0);
    glMatrixMode(GL_MODELVIEW);
    glDisable(GL_DEPTH_TEST);
    glLoadIdentity();
}
/***
 *  *** The "expose_event" signal handler. All the OpenGL re-drawing should
 *   *** be done here. This is repeatedly called as the painting routine
 *    *** every time the 'expose'/'draw' event is signalled.
 *     ***/
static gboolean
expose_event (GtkWidget      *widget,
              GdkEventExpose *,
              gpointer        )
{
    GdkGLContext *glcontext = gtk_widget_get_gl_context (widget);
    GdkGLDrawable *gldrawable = gtk_widget_get_gl_drawable (widget);

    GLfloat width = widget->allocation.width;
    GLfloat height = widget->allocation.height;

  /*** OpenGL BEGIN ***/
    if (!gdk_gl_drawable_gl_begin (gldrawable, glcontext))
	return FALSE;

  if (width != da_width || height != da_height)
	realize(widget,NULL);

    draw_routine (width,height);
  /* Swap buffers */
    gdk_gl_drawable_swap_buffers (gldrawable);

  gdk_gl_drawable_gl_end (gldrawable);
  /*** OpenGL END ***/

  return true;
}

static void
quit_cb (GtkWidget* , gpointer )
{
  gtk_main_quit();
}

void destroy_widget( GtkWidget *, gpointer data)
{
	GtkWidget *widget_to_kill = (GtkWidget *) data;

	gtk_widget_destroy( widget_to_kill );
}

void destroy( void )
{
	if (theState != NULL)
	{
	    theState->pauseFlag = 0;
	    theState->stopFlag = 1;
	}

	gtk_main_quit();
}

static void attr_list_insert( PangoAttrList *attrlist, PangoAttribute *attr )
{
	attr->start_index = 0;
	attr->end_index = G_MAXINT;
	pango_attr_list_insert( attrlist, attr );
}

void select_font(GtkWidget *, gpointer label)
{

  GtkResponseType result;

  GtkWidget *dialog = gtk_font_selection_dialog_new("Select Font");
  do
  {
      result = (GtkResponseType)gtk_dialog_run(GTK_DIALOG(dialog));

      if (result == GTK_RESPONSE_OK || result == GTK_RESPONSE_APPLY)
      {

	PangoFontDescription *font_desc;
	PangoFont *font;
	gchar *fontname = gtk_font_selection_dialog_get_font_name(
				GTK_FONT_SELECTION_DIALOG(dialog));

	font_desc = pango_font_description_from_string(fontname);
	font_list_base = glGenLists (128);

	font = gdk_gl_font_use_pango_font (font_desc, 0, 128, font_list_base);
	if (font == NULL)
	{
	    string err_str ="*** Can't use font '";
	    err_str +=fontname;
	    err_str +="' for OpenGL drawing\n";

	    ShowMessageBox("Warning",err_str);

	    g_print (err_str.c_str());
	    glDeleteLists(font_list_base,128);

	    font_list_base = default_font_list_base;
	    font_desc = default_font_desc;
	    font = gdk_gl_font_use_pango_font (font_desc, 0, 128, font_list_base);
	    result = GTK_RESPONSE_CANCEL;
	}
	else
	{
	    g_print ("*** Loaded font '%s'\n", fontname);
	    glDeleteLists(default_font_list_base,128);
	    default_font_list_base = font_list_base;
	    default_font_desc = font_desc;
	    result = GTK_RESPONSE_CANCEL;
	}

	if (label != NULL)
	    gtk_widget_modify_font(GTK_WIDGET(label), font_desc);

	g_free(fontname);
      }
  }
  while(result == GTK_RESPONSE_OK || result == GTK_RESPONSE_APPLY);

  gtk_widget_destroy(dialog);
}



/**************************************************************************
 *  * The following section contains the GUI building function definitions.
 *   **************************************************************************/

void about_window( void )
{
	GtkWidget *about_window;
	GtkWidget *button;
	GtkWidget *frame;
	GtkWidget *hseparator;
	GtkWidget *label;
	GtkWidget *vbox, *hbox;
	GdkPixmap *apixmap;
	GdkBitmap *mask;
	GtkWidget *aimage;
	PangoAttrList *attrlist;

	about_window = gtk_window_new( GTK_WINDOW_TOPLEVEL );
	gtk_window_set_title( GTK_WINDOW(about_window), "About hiv_sim");
	g_signal_connect( GTK_OBJECT(about_window), "destroy", G_CALLBACK(destroy_widget), about_window );

	vbox = gtk_vbox_new( FALSE, 0 );
	gtk_container_set_border_width( GTK_CONTAINER (vbox), 15);
	gtk_container_add( GTK_CONTAINER( about_window ), vbox );

	button = gtk_button_new_with_label( "OK" );
	g_signal_connect( GTK_OBJECT( button ), "clicked", G_CALLBACK( destroy_widget ), about_window );
	gtk_box_pack_end( GTK_BOX( vbox ), button, FALSE, FALSE, 0 );

	hseparator = gtk_hseparator_new();
	gtk_box_pack_end( GTK_BOX( vbox ), hseparator, TRUE, FALSE, 10 );

	hbox = gtk_hbox_new( FALSE, 10 );
	gtk_box_pack_start( GTK_BOX( vbox ), hbox, FALSE, FALSE, 0 );

	apixmap = gdk_pixmap_create_from_xpm_d( main_window->window, &mask, NULL, hiv_xpm );
	aimage = gtk_image_new_from_pixmap( apixmap, mask );
	gtk_container_add( GTK_CONTAINER( hbox ), aimage );

	vbox = gtk_vbox_new( FALSE, 5 );
	gtk_box_pack_start( GTK_BOX( hbox ), vbox, FALSE, FALSE, 0 );

	attrlist = pango_attr_list_new();
	attr_list_insert( attrlist, pango_attr_size_new( 20*PANGO_SCALE ) );
	attr_list_insert( attrlist, pango_attr_weight_new( PANGO_WEIGHT_BOLD ) );
	attr_list_insert( attrlist, pango_attr_foreground_new( 8738, 39835, 5911 ) );
	label = gtk_label_new( DEFAULT_TITLE " - Version " PACKAGE_VERSION );
	gtk_label_set_attributes( GTK_LABEL(label), attrlist );
	gtk_label_set_justify( GTK_LABEL(label), GTK_JUSTIFY_CENTER );
	gtk_label_set_line_wrap( GTK_LABEL(label), FALSE );
	gtk_box_pack_start( GTK_BOX( vbox ), label, FALSE, FALSE, 0 );
	pango_attr_list_unref( attrlist );

	frame = gtk_frame_new(NULL);
	gtk_frame_set_shadow_type( GTK_FRAME( frame ), GTK_SHADOW_IN );
	gtk_box_pack_start( GTK_BOX( vbox ), frame, FALSE, FALSE, 0 );	

	vbox = gtk_vbox_new( FALSE, 5 );
	gtk_container_set_border_width( GTK_CONTAINER (vbox), 10);
	gtk_container_add( GTK_CONTAINER( frame ), vbox );

	label = gtk_label_new( "A Visualization Program for early HIV infection\nDistributed under the GNU GPL v2.0" );
	gtk_label_set_line_wrap( GTK_LABEL(label), TRUE);
	gtk_box_pack_start( GTK_BOX( vbox ), label, FALSE, FALSE, 0 );	

	label = gtk_label_new( "Authors:\n  David Swan <"PACKAGE_BUGREPORT">\n  Josh Schiffer\n  Dan Reeves" );
	gtk_label_set_justify( GTK_LABEL(label), GTK_JUSTIFY_LEFT);
	gtk_label_set_line_wrap( GTK_LABEL(label), TRUE);
	gtk_box_pack_start( GTK_BOX( vbox ), label, FALSE, FALSE, 0 );

	gtk_widget_show_all(about_window);
}

void sim_ctrl_window( void )
{
	GtkWidget *ctrl_window;
	GtkWidget *button;
	GtkWidget *frame;
	GtkWidget *hseparator;
	GtkWidget *label;
	GtkWidget *vbox, *hbox;
	string ctrl_desc[] = { "Forward: start or resume a simulation run",
  				"Stop: stop a simulation run",
				"Pause: pause a simulation run",
				"Scroll bars: small shifts in time scale when paused or stopped",
				"Zoom in: divide x-axis time range in half",
				"Zoom out: double x-axis time range",
				"Left/Right Arrows: shifts time scale one screen when paused or stopped",
				"Screenshot: save screen as a PNG file",
  				"Save animation: start screen captures at specified interval (for movie creation)"};

	string disp_desc[] = { "Follow top strains: plot virus for strains EVER in the top <n>",
  				"Plot top strains: plot virus for strains currently in the top <n>",
  				"Plot first above cell threshold: plot strains with cells about threshold <n>",
  				"Num strains to plot: used for first two display options",
  				"Log strain cells threshold: used for third display option",
  				"Color by Hamming Distance: strain color based on distance to founder vs number",
  				"Plot Vt vs Time: Include traces of virus by strain and overall",
  				"Plot Actives vs Time: Include traces of actively infected cells",
  				"Plot CD8s vs Time: Include traces of CD8s by strain (if applicable) and overall",
				"Show All Settings: display parameter values (set by Properties tab)",
				"Show Sim Stats: display counts of various outputs (log virus, latent & active cells, etc.)",
				"Enable data writing: required for writing data to files (use -w option to specify the files)",
  				"Scroll time axis: incrementally shift plot as time increases"};

	ctrl_window = gtk_window_new( GTK_WINDOW_TOPLEVEL );
	gtk_window_set_title( GTK_WINDOW(ctrl_window), "Simulation Controls");
	g_signal_connect( GTK_OBJECT(ctrl_window), "destroy", G_CALLBACK(destroy_widget), ctrl_window );

	vbox = gtk_vbox_new( FALSE, 0 );
	gtk_container_set_border_width( GTK_CONTAINER (vbox), 15);
	gtk_container_add( GTK_CONTAINER( ctrl_window ), vbox );

	button = gtk_button_new_with_label( "OK" );
	g_signal_connect( GTK_OBJECT( button ), "clicked", G_CALLBACK( destroy_widget ), ctrl_window );
	gtk_box_pack_end( GTK_BOX( vbox ), button, FALSE, FALSE, 0 );

	hseparator = gtk_hseparator_new();
	gtk_box_pack_end( GTK_BOX( vbox ), hseparator, TRUE, FALSE, 10 );

	hbox = gtk_hbox_new( FALSE, 10 );
	gtk_box_pack_start( GTK_BOX( vbox ), hbox, FALSE, FALSE, 0 );

	vbox = gtk_vbox_new( FALSE, 5 );
	gtk_box_pack_start( GTK_BOX( hbox ), vbox, FALSE, FALSE, 0 );

	frame = gtk_frame_new(NULL);
	gtk_frame_set_shadow_type( GTK_FRAME( frame ), GTK_SHADOW_IN );
	gtk_box_pack_start( GTK_BOX( vbox ), frame, FALSE, FALSE, 0 );	

	vbox = gtk_vbox_new( FALSE, 5 );
	gtk_container_set_border_width( GTK_CONTAINER (vbox), 10);
	gtk_container_add( GTK_CONTAINER( frame ), vbox );

	label = gtk_label_new( "There are several control buttons along the bottom of the display for starting, stopping and pausing the simulation run." );
	gtk_label_set_line_wrap( GTK_LABEL(label), TRUE);
	gtk_box_pack_start( GTK_BOX( vbox ), label, TRUE, FALSE, 0 );	
	hseparator = gtk_hseparator_new();
	gtk_box_pack_start( GTK_BOX( vbox ), hseparator, TRUE, FALSE, 10 );

	for (int i=0; i < sizeof(ctrl_desc)/sizeof(string); i++)
	{
	    label = gtk_label_new(ctrl_desc[i].c_str());
	    gtk_label_set_justify( GTK_LABEL(label), GTK_JUSTIFY_LEFT);
	    gtk_label_set_line_wrap( GTK_LABEL(label), TRUE);
	    gtk_box_pack_start( GTK_BOX( vbox ), label, TRUE, FALSE, 0 );
	}
	hseparator = gtk_hseparator_new();
	gtk_box_pack_start( GTK_BOX( vbox ), hseparator, TRUE, FALSE, 10 );

	label = gtk_label_new( "There are also buttons along the side for selecting display options:");
	gtk_label_set_justify( GTK_LABEL(label), GTK_JUSTIFY_LEFT);
	gtk_label_set_line_wrap( GTK_LABEL(label), TRUE);
	gtk_box_pack_start( GTK_BOX( vbox ), label, TRUE, FALSE, 0 );
	hseparator = gtk_hseparator_new();
	gtk_box_pack_start( GTK_BOX( vbox ), hseparator, TRUE, FALSE, 10 );
	for (int i=0; i < sizeof(disp_desc)/sizeof(string); i++)
	{
	    label = gtk_label_new(disp_desc[i].c_str());
	    gtk_label_set_justify( GTK_LABEL(label), GTK_JUSTIFY_LEFT);
	    gtk_label_set_line_wrap( GTK_LABEL(label), TRUE);
	    gtk_box_pack_start( GTK_BOX( vbox ), label, TRUE, FALSE, 0 );
	}


	gtk_widget_show_all(ctrl_window);
}

void sim_prop_window( void )
{
	GtkWidget *prop_window;
	GtkWidget *button;
	GtkWidget *frame;
	GtkWidget *hseparator;
	GtkWidget *label;
	GtkWidget *vbox, *hbox;
	string var_desc[] = { "Beta mean: mean infectivity for strains",
			      "Beta max: max infectivity for strains",
  				"Beta k: k for weibull on infectivity dist (if selected)",
  				"junk fract: fraction of virus that will lead to defective integrations",
  				"Vmut fract: mutation rate (SNPs)",
  				"eclipse: duration of non-productive phase for infected cells",
  				"Gamma: viral clearance rate",
  				"Pi: viral production per cell per day",
  				"ART eff: ART efficiency (once started)",
  				"ART start: day to start ART",
  				"ART stop: day to suspend ART",
  				"Immun_model: 0=none, 1=global, 2=strain-based, 3=combined",
  				"dImmun: global killing rate per CD8",
  				"dImmun_s: strain-wise killing rate",
  				"theta: scale factor for waxing (infusion) of CD8s",
  				"delta: CD8 death rate",
  				"kappa: scale factor for waning of CD8s",
  				"cd8 mean: used to seed CD8s per strain (model 2)",
  				"cd8 k: used to seed CD8s if weibull",
  				"cd8_0: used to seed CD8s per strain (model 2)",
  				"wax IC50: used to control growth of CD8 population",
  				"wane IC50: used to set carrying capacity of CD8s",
  				"lifespan decay: for increasing immune pressure based on strain life time",
  				"dA IC50: IC50 for Immune model by infected cells",
  				"dA scale: scale factor for increasing immune pressure based on infected cells",
  				"Max Hamming Dist: max groups for strain groups by Hamming Dist",
  				"Delta Hamming Dist: granularity for strain groups by Hamming Dist"};

	prop_window = gtk_window_new( GTK_WINDOW_TOPLEVEL );
	gtk_window_set_title( GTK_WINDOW(prop_window), "Simulation Properties");
	g_signal_connect( GTK_OBJECT(prop_window), "destroy", G_CALLBACK(destroy_widget), prop_window );

	vbox = gtk_vbox_new( FALSE, 0 );
	gtk_container_set_border_width( GTK_CONTAINER (vbox), 15);
	gtk_container_add( GTK_CONTAINER( prop_window ), vbox );

	button = gtk_button_new_with_label( "OK" );
	g_signal_connect( GTK_OBJECT( button ), "clicked", G_CALLBACK( destroy_widget ), prop_window );
	gtk_box_pack_end( GTK_BOX( vbox ), button, FALSE, FALSE, 0 );

	hseparator = gtk_hseparator_new();
	gtk_box_pack_end( GTK_BOX( vbox ), hseparator, TRUE, FALSE, 10 );

	hbox = gtk_hbox_new( FALSE, 10 );
	gtk_box_pack_start( GTK_BOX( vbox ), hbox, FALSE, FALSE, 0 );

	vbox = gtk_vbox_new( FALSE, 5 );
	gtk_box_pack_start( GTK_BOX( hbox ), vbox, FALSE, FALSE, 0 );

	frame = gtk_frame_new(NULL);
	gtk_frame_set_shadow_type( GTK_FRAME( frame ), GTK_SHADOW_IN );
	gtk_box_pack_start( GTK_BOX( vbox ), frame, FALSE, FALSE, 0 );	

	vbox = gtk_vbox_new( FALSE, 5 );
	gtk_container_set_border_width( GTK_CONTAINER (vbox), 10);
	gtk_container_add( GTK_CONTAINER( frame ), vbox );

	label = gtk_label_new( "There are text entry boxes for many of the simulation parameters.");
	gtk_label_set_justify( GTK_LABEL(label), GTK_JUSTIFY_LEFT);
	gtk_box_pack_start( GTK_BOX( vbox ), label, TRUE, FALSE, 0 );
	label = gtk_label_new( "Changes to these take effect when a new run is started or the 'Apply' button is pressed." );
	gtk_label_set_justify( GTK_LABEL(label), GTK_JUSTIFY_LEFT);
	gtk_box_pack_start( GTK_BOX( vbox ), label, TRUE, FALSE, 0 );

	hseparator = gtk_hseparator_new();
	gtk_box_pack_start( GTK_BOX( vbox ), hseparator, TRUE, FALSE, 10 );

	for (int i=0; i < sizeof(var_desc)/sizeof(string); i++)
	{
	    label = gtk_label_new(var_desc[i].c_str());
	    gtk_label_set_justify( GTK_LABEL(label), GTK_JUSTIFY_LEFT);
	    gtk_box_pack_start( GTK_BOX( vbox ), label, TRUE, TRUE, 0 );
	}

	hseparator = gtk_hseparator_new();
	gtk_box_pack_start( GTK_BOX( vbox ), hseparator, TRUE, FALSE, 10 );

	label = gtk_label_new( "There are also options for changing the background color, default font and max range for y-axis values.");
	gtk_label_set_justify( GTK_LABEL(label), GTK_JUSTIFY_LEFT);
	gtk_label_set_line_wrap( GTK_LABEL(label), TRUE);
	gtk_box_pack_start( GTK_BOX( vbox ), label, TRUE, FALSE, 0 );

	gtk_widget_show_all(prop_window);
}

void make_gui(void)
{
#ifdef __sun
	static GtkActionEntry menu_entries[] = {
		{ "FileMenu", NULL, "_File" },
		{ "HelpMenu", NULL, "_Help" },
		{ "Quit", GTK_STOCK_QUIT, "_Quit", "<control>Q", "Exit the program", destroy },
		{ "HelpAbout", NULL, "_About", NULL, "About hiv_sim...", about_window },
		{ "HelpSim", NULL, "_Simulation", NULL, "Simulation Controls...", sim_ctrl_window },
		{ "HelpProperties", NULL, "_Properties", NULL, "Simulation Properties...", sim_prop_window },
	};
#else
	static GtkActionEntry menu_entries[] = {
		{ "FileMenu", NULL, "_File" },
		{ "HelpMenu", NULL, "_Help" },
		{ "Quit", GTK_STOCK_QUIT, "_Quit", "<control>Q", "Exit the program", destroy },
		{ "HelpAbout", GTK_STOCK_ABOUT, "_About", NULL, "About hiv_sim...", about_window },
		{ "HelpSim", NULL, "_Simulation", NULL, "Simulation Controls...", sim_ctrl_window },
		{ "HelpProperties", NULL, "_Properties", NULL, "Simulation Properties...", sim_prop_window },
	};
#endif
	static const char *ui_description =
	"<ui>"
	"  <menubar name='MainMenu'>"
	"    <menu action='FileMenu'>"
	"      <menuitem action='Quit'/>"
	"    </menu>"
	"    <menu action='HelpMenu'>"
	"      <menuitem action='HelpAbout'/>"
	"      <menuitem action='HelpSim'/>"
	"      <menuitem action='HelpProperties'/>"
	"    </menu>"
	"  </menubar>"
	"</ui>";
	GtkActionGroup *action_group;
	GtkUIManager *ui_manager;
	GtkAccelGroup *accel_group;
	GError *error;
	GtkWidget *menubar;
	GtkWidget *moviebutton = NULL;

	GtkObject *spinbutton_adj;
	GtkObject *strain_spinbutton_adj;
	GtkObject *cell_spinbutton_adj;

	GtkWidget *button;
	GtkWidget *frame;
	GtkWidget *label;
	GtkWidget *notebook;
	GtkWidget *vbox,*svbox,*hbox;
	GtkWidget *optBox;

	GtkWidget *writeEnable;
	GtkWidget *scrollAxes;

	GtkWidget *plot_vt_button;
	GtkWidget *plot_act_button;
	GtkWidget *plot_CD8_button;
	GtkWidget *show_params_button;
	GtkWidget *show_stats_button;

	/* Creates the main Window ; sets its title ; attaches the closing behavior */
	main_window = gtk_window_new( GTK_WINDOW_TOPLEVEL );
	gtk_window_set_default_size( GTK_WINDOW(main_window), DEFAULT_WIDTH,DEFAULT_HEIGHT);
	gtk_window_set_title( GTK_WINDOW(main_window), DEFAULT_TITLE);
	g_signal_connect( GTK_OBJECT(main_window), "destroy", G_CALLBACK(destroy), NULL );
	gtk_quit_add_destroy(1, GTK_OBJECT(main_window));
        /* Connect signal handlers to the window */
        g_signal_connect (G_OBJECT (main_window), "delete_event", G_CALLBACK (quit_cb), NULL);
	vbox = gtk_vbox_new( FALSE, 0 );
	gtk_container_add( GTK_CONTAINER( main_window ), vbox );

	/* Create Menu */
	action_group = gtk_action_group_new ("MenuActions");
	gtk_action_group_add_actions( action_group, menu_entries, G_N_ELEMENTS(menu_entries), main_window);
	ui_manager = gtk_ui_manager_new();
	gtk_ui_manager_insert_action_group( ui_manager, action_group, 0 );
	accel_group = gtk_ui_manager_get_accel_group( ui_manager );
	gtk_window_add_accel_group( GTK_WINDOW(main_window), accel_group );
	error = NULL;
	if( !gtk_ui_manager_add_ui_from_string( ui_manager, ui_description, -1, &error ) ) {
		g_message( "Building menus failed: %s", error->message );
		g_error_free( error );
		exit( EXIT_FAILURE );
	}
	menubar = gtk_ui_manager_get_widget( ui_manager, "/MainMenu" );
	gtk_box_pack_start( GTK_BOX( vbox ), menubar, FALSE, FALSE, 0 );

	/* Create the Notebook */
	notebook = gtk_notebook_new();
	gtk_notebook_set_tab_pos( GTK_NOTEBOOK( notebook ), GTK_POS_TOP );
	gtk_box_pack_start( GTK_BOX( vbox ), notebook, TRUE, TRUE, 0 );


	/* Create Simulation TAB */
	vbox = gtk_vbox_new( FALSE, 10 );
	gtk_container_set_border_width( GTK_CONTAINER( vbox ), 10 );
	gtk_notebook_append_page( GTK_NOTEBOOK( notebook ), vbox, gtk_label_new( "Simulation" ) );
	/* Create frame to hold the graph */
	hbox = gtk_hbox_new( FALSE, 0 );
	gtk_box_pack_start( GTK_BOX( vbox ), hbox, TRUE, TRUE, 0 );
	//frame = gtk_frame_new(NULL);
	//gtk_frame_set_shadow_type( GTK_FRAME(frame), GTK_SHADOW_IN );
	//gtk_box_pack_start( GTK_BOX( hbox ), frame, TRUE, TRUE, 0 );
	gtk_widget_show (hbox);

	/* Create glconfig attributes for gtkglext */
	glconfig = gdk_gl_config_new_by_mode( (GdkGLConfigMode)(GDK_GL_MODE_RGB | GDK_GL_MODE_DOUBLE ));
	if( glconfig == NULL ) {
		g_print ("*** Problem with gtkglext.\n");
		exit(1);
	}
	image = gtk_drawing_area_new();
	gtk_widget_set_size_request( GTK_WIDGET(image), DEFAULT_WIDTH, DEFAULT_HEIGHT);
	gtk_widget_set_gl_capability( image, glconfig, NULL, TRUE, GDK_GL_RGBA_TYPE );
	/* Connect signal handlers to the drawing area */
	g_signal_connect_after (GTK_OBJECT (image), "realize",
				G_CALLBACK (realize), NULL);
	g_signal_connect_after (GTK_OBJECT (image), "configure_event",
				G_CALLBACK (configure_event), NULL);
	
	g_signal_connect(GTK_OBJECT(image), "expose_event", G_CALLBACK(expose_event), NULL);

	gtk_box_pack_start( GTK_BOX( hbox ), image, TRUE, TRUE, 0 );
	//gtk_container_add( GTK_CONTAINER(frame), image );

	/* Create vbox to hold plot settings */
	optBox = gtk_vbox_new( FALSE, 0 );
	gtk_box_pack_end( GTK_BOX(hbox), optBox, FALSE, FALSE, 0 );
	gtk_widget_show (optBox);


	followTopStrains = gtk_check_button_new_with_label( "Follow top strains" );
	g_signal_connect( GTK_OBJECT( followTopStrains ), "toggled", G_CALLBACK(setCheckFollowStrains), NULL);
	gtk_box_pack_start( GTK_BOX( optBox ), followTopStrains, TRUE, FALSE, 0 );

	plotTopStrains = gtk_check_button_new_with_label( "Plot top strains" );
	g_signal_connect( GTK_OBJECT( plotTopStrains ), "toggled", G_CALLBACK(setCheckTopStrains), NULL);
	gtk_box_pack_start( GTK_BOX( optBox ), plotTopStrains, TRUE, FALSE, 0 );

	plotFirstStrains = gtk_check_button_new_with_label( "Plot first above cell threshold" );
	g_signal_connect( GTK_OBJECT( plotFirstStrains ), "toggled", G_CALLBACK(setCheckFirstStrains), NULL);
	gtk_box_pack_start( GTK_BOX( optBox ), plotFirstStrains, TRUE, FALSE, 0 );

	top_strain_box = gtk_hbox_new( FALSE, 0 );
	gtk_box_pack_start( GTK_BOX(optBox), top_strain_box, TRUE, FALSE, 0 );
	gtk_widget_show (top_strain_box);

	label = gtk_label_new( "Num strains to plot:" );
	gtk_box_pack_start( GTK_BOX( top_strain_box ), label, FALSE, FALSE, 0);
	strain_spinbutton_adj = gtk_adjustment_new( 1, 1, MAX_STRAIN_COLORS, 1, 1, 0 );

	strain_widget = gtk_spin_button_new( GTK_ADJUSTMENT(strain_spinbutton_adj), 1, 0 );
	gtk_box_pack_start( GTK_BOX( top_strain_box ), strain_widget, FALSE, FALSE, 0);
	g_signal_connect (G_OBJECT (strain_widget), "value_changed",
			  G_CALLBACK (setMaxStrains), NULL);

	first_strain_box = gtk_hbox_new( FALSE, 0 );
	gtk_box_pack_start( GTK_BOX(optBox), first_strain_box, TRUE, FALSE, 0 );
	gtk_widget_show (first_strain_box);

	label = gtk_label_new( "Log strain cells threshold:" );
	gtk_box_pack_start( GTK_BOX( first_strain_box ), label, FALSE, FALSE, 0);
	cell_spinbutton_adj = gtk_adjustment_new( 1, 1, 5, 1, 1, 0 );

	cell_widget = gtk_spin_button_new( GTK_ADJUSTMENT(cell_spinbutton_adj), 1, 0 );
	gtk_box_pack_start( GTK_BOX( first_strain_box ), cell_widget, FALSE, FALSE, 0);
	g_signal_connect (G_OBJECT (cell_widget), "value_changed",
			  G_CALLBACK (setCellThreshold), NULL);


	colorByHamming = gtk_check_button_new_with_label( "Color by Hamming Distance" );
	g_signal_connect( GTK_OBJECT( colorByHamming ), "toggled", G_CALLBACK(setCheckColorByHamming), NULL);
	gtk_box_pack_start( GTK_BOX( optBox ), colorByHamming, TRUE, FALSE, 0 );

	plot_vt_button = gtk_check_button_new_with_label( "Plot Vt vs time" );
	g_signal_connect( GTK_OBJECT( plot_vt_button ), "toggled", G_CALLBACK(setCheckVt), NULL);
	gtk_box_pack_start( GTK_BOX( optBox ), plot_vt_button, TRUE, FALSE, 0 );

	plot_act_button = gtk_check_button_new_with_label( "Plot Actives vs time" );
	g_signal_connect( GTK_OBJECT( plot_act_button ), "toggled", G_CALLBACK(setCheckAct), NULL);
	gtk_box_pack_start( GTK_BOX( optBox ), plot_act_button, TRUE, FALSE, 0 );

	plot_CD8_button = gtk_check_button_new_with_label( "Plot CD8s vs time" );
	g_signal_connect( GTK_OBJECT( plot_CD8_button ), "toggled", G_CALLBACK(setCheckCD8s), NULL);
	gtk_box_pack_start( GTK_BOX( optBox ), plot_CD8_button, TRUE, FALSE, 0 );

	show_params_button = gtk_check_button_new_with_label( "Show All Settings" );
	g_signal_connect( GTK_OBJECT( show_params_button ), "toggled", G_CALLBACK(setShowParams), NULL);
	gtk_box_pack_start( GTK_BOX( optBox ), show_params_button, TRUE, FALSE, 0 );

	show_stats_button = gtk_check_button_new_with_label( "Show Sim Stats" );
	g_signal_connect( GTK_OBJECT( show_stats_button ), "toggled", G_CALLBACK(setShowStats), NULL);
	gtk_box_pack_start( GTK_BOX( optBox ), show_stats_button, TRUE, FALSE, 0 );

	writeEnable = gtk_check_button_new_with_label( "Enable data writing" );
	g_signal_connect( GTK_OBJECT( writeEnable ), "toggled", G_CALLBACK(setCheckWriteOn), NULL);
	gtk_box_pack_start( GTK_BOX( optBox ), writeEnable, TRUE, FALSE, 0 );

	scrollAxes = gtk_check_button_new_with_label( "Scroll time axes" );
	g_signal_connect( GTK_OBJECT( scrollAxes ), "toggled", G_CALLBACK(setCheckScrollAxes), NULL);
	gtk_box_pack_start( GTK_BOX( optBox ), scrollAxes, TRUE, FALSE, 0 );

	/* Create the bottom box for pause, play, etc. buttons */
	hbox = gtk_hbox_new( FALSE, 10 );
	/*
       *    * Run simulation button.
       *       */

	run_button = gtk_button_new_from_stock( GTK_STOCK_GO_FORWARD );

	g_signal_connect (G_OBJECT (run_button), "clicked",
			  G_CALLBACK (run_cb), NULL);

	gtk_box_pack_start (GTK_BOX (hbox), run_button, TRUE, TRUE, 5);
	
	gtk_widget_show (run_button);
	/*
       *    * Stop simulation button.
       *       */

	stop_button = gtk_button_new_from_stock( GTK_STOCK_STOP );

	g_signal_connect (G_OBJECT (stop_button), "clicked",
			  G_CALLBACK (stop_cb), NULL);

	gtk_box_pack_start (GTK_BOX (hbox), stop_button, TRUE, TRUE, 5);
	
	gtk_widget_show (stop_button);
	gtk_widget_set_sensitive( GTK_WIDGET(stop_button), FALSE );
	/*
       *    * Pause simulation button.
       *       */

	pause_button = gtk_button_new_with_label( "Pause" );

	g_signal_connect (G_OBJECT (pause_button), "clicked",
			  G_CALLBACK (pause_cb), NULL);

	gtk_box_pack_start (GTK_BOX (hbox), pause_button, TRUE, TRUE, 5);
	
	gtk_widget_show (pause_button);
	gtk_widget_set_sensitive( GTK_WIDGET(pause_button), FALSE );

	/* Create cell time scroll bar (hidden for now) */
        gdouble lower = 0.0;
        gdouble upper = 20.01;
        gdouble step_increment = 0.01;
        gdouble page_increment = 0.1;
        gdouble page_size = 0.01;

	adjuster = gtk_adjustment_new (20.0, lower, upper, step_increment,
                                                         page_increment, page_size);


	backup_slider = gtk_hscrollbar_new(GTK_ADJUSTMENT(adjuster));

	gtk_box_pack_start (GTK_BOX (hbox), backup_slider, TRUE, TRUE, 5);

	g_signal_connect (G_OBJECT (adjuster), "value-changed",
			  G_CALLBACK (scrollTime), NULL);
	/*
       *    * Zoom in simulation button.
       *       */

	zoomin_button = gtk_button_new_with_label ("Zoom in");

	g_signal_connect (G_OBJECT (zoomin_button), "clicked",
			  G_CALLBACK (zoomin_cb), NULL);

	gtk_box_pack_start (GTK_BOX (hbox), zoomin_button, TRUE, TRUE, 5);
	
	gtk_widget_show (zoomin_button);

	/*
       *    * Zoom out simulation button.
       *       */

	zoomout_button = gtk_button_new_with_label ("Zoom out");

	g_signal_connect (G_OBJECT (zoomout_button), "clicked",
			  G_CALLBACK (zoomout_cb), NULL);

	gtk_box_pack_start (GTK_BOX (hbox), zoomout_button, TRUE, TRUE, 5);
	
	gtk_widget_show (zoomout_button);

	gtk_widget_set_sensitive( GTK_WIDGET(zoomin_button), TRUE );
	gtk_widget_set_sensitive( GTK_WIDGET(zoomout_button), TRUE );
	/*
       *    * Page left simulation button.
       *       */

	page_left_button = gtk_button_new_with_label ("<<");

	g_signal_connect (G_OBJECT (page_left_button), "clicked",
			  G_CALLBACK (page_left_cb), NULL);

	gtk_box_pack_start (GTK_BOX (hbox), page_left_button, TRUE, TRUE, 5);
	
	gtk_widget_show (page_left_button);

	/*
       *    * Page right simulation button.
       *       */

	page_right_button = gtk_button_new_with_label (">>");

	g_signal_connect (G_OBJECT (page_right_button), "clicked",
			  G_CALLBACK (page_right_cb), NULL);

	gtk_box_pack_start (GTK_BOX (hbox), page_right_button, TRUE, TRUE, 5);
	
	gtk_widget_show (page_right_button);
	gtk_widget_set_sensitive( GTK_WIDGET(page_left_button), FALSE );
	gtk_widget_set_sensitive( GTK_WIDGET(page_right_button), FALSE );

	/* Button to take screenshot */
	button = gtk_button_new_with_label( "Screenshot" );
	g_signal_connect( GTK_OBJECT( button ), "clicked",G_CALLBACK( open_file_selection ), (void *)&take_screenshot );
	gtk_box_pack_start( GTK_BOX( hbox ), button, TRUE, TRUE, 5 );
	gtk_widget_show (button);
	/* Record to movie */
	moviebutton = gtk_check_button_new_with_label( "Save animation." );
	g_signal_connect( GTK_OBJECT( moviebutton ), "clicked", G_CALLBACK(check_button_clicked), (void *)("movie"));
	//gtk_widget_set_sensitive( GTK_WIDGET(moviebutton), FALSE );
	gtk_box_pack_start( GTK_BOX( hbox ), moviebutton, TRUE, TRUE, 5 );
	gtk_widget_show (moviebutton);
	gtk_widget_show (hbox);
	gtk_box_pack_end( GTK_BOX( vbox ), hbox, FALSE, FALSE, 0 );

	/* Create Properties TAB */
	vbox = gtk_vbox_new( FALSE, 0 );
	gtk_notebook_append_page( GTK_NOTEBOOK( notebook ), vbox, gtk_label_new( "Properties" ) );
	gtk_container_set_border_width( GTK_CONTAINER( vbox ), 10 );

	/* Parameters frame */
	frame = gtk_frame_new( "Simulation Parameters" );
	gtk_box_pack_start( GTK_BOX( vbox ), frame, FALSE, FALSE, 0);

  	GtkWidget *labels[MAX_GUI_ENTRIES];
	int label_num=0;
  	labels[label_num++]= gtk_label_new("Beta mean:");
  	labels[label_num++]= gtk_label_new("Beta max:");
  	labels[label_num++]= gtk_label_new("Beta k:");

  	labels[label_num++]= gtk_label_new("junk fract:");
  	labels[label_num++]= gtk_label_new("Vmut fract:");
  	labels[label_num++]= gtk_label_new("eclipse:");

  	labels[label_num++]= gtk_label_new("Gam:");
  	labels[label_num++]= gtk_label_new("Pi:");

  	labels[label_num++]= gtk_label_new("ART eff:");
  	labels[label_num++]= gtk_label_new("ART start:");
  	labels[label_num++]= gtk_label_new("ART stop:");

  	labels[label_num++]= gtk_label_new("Immun_model:");
  	labels[label_num++]= gtk_label_new("dImmun:");
  	labels[label_num++]= gtk_label_new("dImmun_s:");

  	labels[label_num++]= gtk_label_new("theta:");
  	labels[label_num++]= gtk_label_new("delta:");
  	labels[label_num++]= gtk_label_new("kappa:");

  	labels[label_num++]= gtk_label_new("cd8 mean:");
  	labels[label_num++]= gtk_label_new("cd8 k:");
  	labels[label_num++]= gtk_label_new("cd8_0:");

  	labels[label_num++]= gtk_label_new("wax ic50:");
  	labels[label_num++]= gtk_label_new("wane ic50:");
  	labels[label_num++]= gtk_label_new("lifespan decay:");

  	labels[label_num++]= gtk_label_new("dA_ic50:");
  	labels[label_num++]= gtk_label_new("dA scale:");

  	labels[label_num++]= gtk_label_new("Max Hamming Dist:");
  	labels[label_num++]= gtk_label_new("Delta Hamming Dist:");

  	GtkWidget *entry[MAX_GUI_ENTRIES];
	int entry_num=0;

	beta_entry = gtk_entry_new();
	entry[entry_num++]=beta_entry;
	beta_max_entry = gtk_entry_new();
	entry[entry_num++]=beta_max_entry;
	beta_k_entry = gtk_entry_new();
	entry[entry_num++]=beta_k_entry;

	junk_entry = gtk_entry_new();
	entry[entry_num++]=junk_entry;
	mu_entry = gtk_entry_new();
	entry[entry_num++]=mu_entry;
	eclipse_entry = gtk_entry_new();
	entry[entry_num++]=eclipse_entry;

	gam_entry = gtk_entry_new();
	entry[entry_num++]=gam_entry;
	pi_entry = gtk_entry_new();
	entry[entry_num++]=pi_entry;

	ART_eff_entry = gtk_entry_new();
	entry[entry_num++]=ART_eff_entry;
	ART_start_entry = gtk_entry_new();
	entry[entry_num++]=ART_start_entry;
	ART_stop_entry = gtk_entry_new();
	entry[entry_num++]=ART_stop_entry;

	Immun_model_entry = gtk_entry_new();
	entry[entry_num++]=Immun_model_entry;
	dImmun_entry = gtk_entry_new();
	entry[entry_num++]=dImmun_entry;
	dImmun_s_entry = gtk_entry_new();
	entry[entry_num++]=dImmun_s_entry;

	theta_entry = gtk_entry_new();
	entry[entry_num++]=theta_entry;
	delta_entry = gtk_entry_new();
	entry[entry_num++]=delta_entry;
	kappa_entry = gtk_entry_new();
	entry[entry_num++]=kappa_entry;

	cd8_mean_entry = gtk_entry_new();
	entry[entry_num++]=cd8_mean_entry;
	cd8_0_entry = gtk_entry_new();
	entry[entry_num++]=cd8_0_entry;
	cd8_k_entry = gtk_entry_new();
	entry[entry_num++]=cd8_k_entry;

	dImmun_IC50_entry = gtk_entry_new();
	entry[entry_num++]=dImmun_IC50_entry;
	cd8_ic50_entry = gtk_entry_new();
	entry[entry_num++]=cd8_ic50_entry;
	lifespan_decay_entry = gtk_entry_new();
	entry[entry_num++]=lifespan_decay_entry;

	dA_ic50_entry = gtk_entry_new();
	entry[entry_num++]=dA_ic50_entry;
	dA_scale_entry = gtk_entry_new();
	entry[entry_num++]=dA_scale_entry;

	maxHamming_entry = gtk_entry_new();
	entry[entry_num++]=maxHamming_entry;
	deltaHamming_entry = gtk_entry_new();
	entry[entry_num++]=deltaHamming_entry;

	if (label_num != entry_num) {
	    fprintf(stderr,"Internal error %d labels and %d widgets\nExitings!\n",label_num,entry_num);
	    exit(1);
	}
	if (label_num > MAX_GUI_ENTRIES) {
	    fprintf(stderr,"Internal error too many GUI variables (%d vs %d)\nExiting!\n",
		label_num,MAX_GUI_ENTRIES);
	    exit(1);
	}
	GtkWidget *table = gtk_table_new(round(label_num/3.0), 6, FALSE);
        gtk_container_add(GTK_CONTAINER(frame), table);
	gtk_widget_show (table);

	int rownum=0;
	int colnum=0;
	for (int i=0; i < label_num; i++) {
	    gtk_table_attach_defaults(GTK_TABLE(table), labels[i], colnum, colnum+1, rownum, rownum+1);
	    gtk_table_attach_defaults(GTK_TABLE(table), entry[i], colnum+1, colnum+2, rownum, rownum+1);
	    colnum+=2;
	    if ((i+1) % 3 == 0) {
		rownum++;
		colnum=0;
	    }
	}

	char value_str[100];
	sprintf(value_str,"%lf", theState->ART_eff);
	gtk_entry_set_text (GTK_ENTRY(ART_eff_entry), value_str);
	sprintf(value_str,"%lf", theState->ART_start);
	gtk_entry_set_text (GTK_ENTRY(ART_start_entry), value_str);
	sprintf(value_str,"%lf", theState->ART_stop);
	gtk_entry_set_text (GTK_ENTRY(ART_stop_entry), value_str);

	sprintf(value_str,"%e", theState->Bt_mean);
	gtk_entry_set_text (GTK_ENTRY(beta_entry), value_str);
	sprintf(value_str,"%e", theState->Bt_max);
	gtk_entry_set_text (GTK_ENTRY(beta_max_entry), value_str);
	sprintf(value_str,"%lf", theState->Bt_k);
	gtk_entry_set_text (GTK_ENTRY(beta_k_entry), value_str);
	sprintf(value_str,"%lf", theState->lifespan_decay);
	gtk_entry_set_text (GTK_ENTRY(lifespan_decay_entry), value_str);

	sprintf(value_str,"%lf", theState->gam);
	gtk_entry_set_text (GTK_ENTRY(gam_entry), value_str);
	sprintf(value_str,"%lf", theState->mu0);
	gtk_entry_set_text (GTK_ENTRY(mu_entry), value_str);
	sprintf(value_str,"%lf", theState->pi);
	gtk_entry_set_text (GTK_ENTRY(pi_entry), value_str);

	sprintf(value_str,"%lf", theState->theta);
	gtk_entry_set_text (GTK_ENTRY(theta_entry), value_str);
	sprintf(value_str,"%lf", theState->junk);
	gtk_entry_set_text (GTK_ENTRY(junk_entry), value_str);
	sprintf(value_str,"%lf", theState->eclipse);
	gtk_entry_set_text (GTK_ENTRY(eclipse_entry), value_str);
	sprintf(value_str,"%lf", theState->dA_ic50_0);
	gtk_entry_set_text (GTK_ENTRY(dA_ic50_entry), value_str);
	sprintf(value_str,"%lf", theState->dA_scale);
	gtk_entry_set_text (GTK_ENTRY(dA_scale_entry), value_str);
	sprintf(value_str,"%lf", theState->delta);
	gtk_entry_set_text (GTK_ENTRY(delta_entry), value_str);
	sprintf(value_str,"%d", theState->maxTopStrains);
	gtk_entry_set_text (GTK_ENTRY(strain_widget), value_str);
	sprintf(value_str,"%d", theState->cellThreshold);
	gtk_entry_set_text (GTK_ENTRY(cell_widget), value_str);
	sprintf(value_str,"%lf", theState->cd8_mean);
	gtk_entry_set_text (GTK_ENTRY(cd8_mean_entry), value_str);
	sprintf(value_str,"%lf", theState->cd8_0);
	gtk_entry_set_text (GTK_ENTRY(cd8_0_entry), value_str);
	sprintf(value_str,"%lf", theState->cd8_k);
	gtk_entry_set_text (GTK_ENTRY(cd8_k_entry), value_str);
	sprintf(value_str,"%lf", theState->dImmun);
	gtk_entry_set_text (GTK_ENTRY(dImmun_entry), value_str);
	sprintf(value_str,"%lf", theState->dImmun_IC50_0);
	gtk_entry_set_text (GTK_ENTRY(dImmun_IC50_entry), value_str);
	sprintf(value_str,"%lf", theState->dImmun_s);
	gtk_entry_set_text (GTK_ENTRY(dImmun_s_entry), value_str);
	sprintf(value_str,"%lf", theState->kappa);
	gtk_entry_set_text (GTK_ENTRY(kappa_entry), value_str);
	sprintf(value_str,"%lf", theState->cd8_ic50_0);
	gtk_entry_set_text (GTK_ENTRY(cd8_ic50_entry), value_str);
	sprintf(value_str,"%d", theState->immun_model);
	gtk_entry_set_text (GTK_ENTRY(Immun_model_entry), value_str);
	sprintf(value_str,"%d", theState->maxHammingDist);
	gtk_entry_set_text (GTK_ENTRY(maxHamming_entry), value_str);
	sprintf(value_str,"%d", theState->cd8GroupSize);
	gtk_entry_set_text (GTK_ENTRY(deltaHamming_entry), value_str);
	button = gtk_button_new_with_label( "Apply (now vs. next run)..." );
	gtk_box_pack_start( GTK_BOX( vbox ), button, FALSE, FALSE, 0 );
	g_signal_connect( G_OBJECT( button ), "clicked",G_CALLBACK(updateParams), (void *)button);


	/* GUI Option frame */
	frame = gtk_frame_new( "GUI Options" );
	gtk_box_pack_start( GTK_BOX( vbox ), frame, FALSE, FALSE, 0);
	svbox = gtk_vbox_new( FALSE, 10 );
	gtk_container_set_border_width( GTK_CONTAINER( svbox ), 5 );
	gtk_container_add( GTK_CONTAINER( frame ), svbox );

	/* Font name */
	hbox = gtk_hbox_new( FALSE, 10 );
	gtk_box_pack_start( GTK_BOX( svbox ), hbox, FALSE, FALSE, 0);
	label = gtk_label_new( "Font:" );
	gtk_box_pack_start( GTK_BOX( hbox ), label, FALSE, FALSE, 0);
	font_name_widget = gtk_entry_new();
	gtk_box_pack_start( GTK_BOX( hbox ), font_name_widget, TRUE, TRUE, 0 );
	gtk_entry_set_text (GTK_ENTRY(font_name_widget), "Courier bold");

	/* Font size */
	label = gtk_label_new( "Font size:" );
	gtk_box_pack_start( GTK_BOX( hbox ), label, FALSE, FALSE, 0);
	spinbutton_adj = gtk_adjustment_new( 12, 8, 24, 2, 2, 0 );



	font_size_widget = gtk_spin_button_new( GTK_ADJUSTMENT(spinbutton_adj), 1, 0 );
	gtk_box_pack_start( GTK_BOX( hbox ), font_size_widget, FALSE, FALSE, 0);

	button = gtk_button_new_with_label( "Load..." );
	gtk_box_pack_start( GTK_BOX( hbox ), button, FALSE, TRUE, 0 );
	g_signal_connect( G_OBJECT( button ), "clicked",G_CALLBACK(font_changed), (void *)button);

	GtkWidget *font_dialog=gtk_font_selection_dialog_new ("Choose a new font...");
	button = gtk_button_new_with_label( "Browse fonts..." );
	gtk_box_pack_start( GTK_BOX( hbox ), button, FALSE, TRUE, 0 );
	g_signal_connect( G_OBJECT( button ), "clicked",G_CALLBACK(select_font), (void *)font_dialog);

	/* Background RGB */
	hbox = gtk_hbox_new( FALSE, 10 );
	gtk_box_pack_start( GTK_BOX( svbox ), hbox, FALSE, FALSE, 0);
	rgb_box = gtk_event_box_new( );
	label = gtk_label_new( "Background RGB:" );
	gtk_container_add (GTK_CONTAINER (rgb_box), label);
	gtk_box_pack_start( GTK_BOX( hbox ), rgb_box, FALSE, FALSE, 0);

	spinbutton_adj = gtk_adjustment_new( 255, 0, MAX_RGB_VAL, 1, 10, 0 );
	background_r = gtk_spin_button_new( GTK_ADJUSTMENT(spinbutton_adj), 1, 0 );
	gtk_box_pack_start( GTK_BOX( hbox ), background_r, FALSE, FALSE, 0);
	label = gtk_label_new( "R" );
	gtk_box_pack_start( GTK_BOX( hbox ), label, FALSE, FALSE, 0);
	g_signal_connect( G_OBJECT(background_r),"changed",G_CALLBACK(background_changed), (void *)NULL );

	spinbutton_adj = gtk_adjustment_new( 255, 0, MAX_RGB_VAL, 1, 10, 0 );
	background_g = gtk_spin_button_new( GTK_ADJUSTMENT(spinbutton_adj), 1, 0 );
	gtk_box_pack_start( GTK_BOX( hbox ), background_g, FALSE, FALSE, 0);
	label = gtk_label_new( "G" );
	gtk_box_pack_start( GTK_BOX( hbox ), label, FALSE, FALSE, 0);
	g_signal_connect( G_OBJECT(background_g),"changed",G_CALLBACK(background_changed), (void *)NULL );

	spinbutton_adj = gtk_adjustment_new( 255, 0, MAX_RGB_VAL, 1, 10, 0 );
	background_b = gtk_spin_button_new( GTK_ADJUSTMENT(spinbutton_adj), 1, 0 );
	gtk_box_pack_start( GTK_BOX( hbox ), background_b, FALSE, FALSE, 0);
	label = gtk_label_new( "B" );
	gtk_box_pack_start( GTK_BOX( hbox ), label, FALSE, FALSE, 0);
	g_signal_connect( G_OBJECT(background_b),"changed",G_CALLBACK(background_changed), (void *)NULL );

	GdkColor init_color;

	init_color.red = 65535;
	init_color.green = 65535;
	init_color.blue = 65535;

	gtk_widget_modify_bg(rgb_box,GTK_STATE_NORMAL, &init_color);
	gtk_widget_set_sensitive( GTK_WIDGET(backup_slider), TRUE );
	gtk_widget_hide( GTK_WIDGET(backup_slider) );

	hbox = gtk_hbox_new( FALSE, 10 );
	gtk_box_pack_start( GTK_BOX( svbox ), hbox, FALSE, FALSE, 0);
	label = gtk_label_new( "Max log virus" );
	gtk_box_pack_start( GTK_BOX( hbox ), label, FALSE, FALSE, 0);
	max_log_v_entry = gtk_entry_new();
	sprintf(value_str,"%g", theState->max_log_v_value);
	gtk_entry_set_text (GTK_ENTRY(max_log_v_entry), value_str);
	gtk_box_pack_start( GTK_BOX( hbox ), max_log_v_entry, FALSE, FALSE, 0);

	hbox = gtk_hbox_new( FALSE, 10 );
	gtk_box_pack_start( GTK_BOX( svbox ), hbox, FALSE, FALSE, 0);
	label = gtk_label_new( "Max log A cells" );
	gtk_box_pack_start( GTK_BOX( hbox ), label, FALSE, FALSE, 0);
	max_log_a_entry = gtk_entry_new();
	sprintf(value_str,"%g", theState->max_log_a_value);
	gtk_entry_set_text (GTK_ENTRY(max_log_a_entry), value_str);
	gtk_box_pack_start( GTK_BOX( hbox ), max_log_a_entry, FALSE, FALSE, 0);

	hbox = gtk_hbox_new( FALSE, 10 );
	gtk_box_pack_start( GTK_BOX( svbox ), hbox, FALSE, FALSE, 0);
	label = gtk_label_new( "Max cd8 cells" );
	gtk_box_pack_start( GTK_BOX( hbox ), label, FALSE, FALSE, 0);
	max_log_cd8_entry = gtk_entry_new();
	sprintf(value_str,"%lf", theState->max_log_cd8_value);
	gtk_entry_set_text (GTK_ENTRY(max_log_cd8_entry), value_str);
	gtk_box_pack_start( GTK_BOX( hbox ), max_log_cd8_entry, FALSE, FALSE, 0);

	gtk_widget_show (max_log_v_entry);
	gtk_widget_show (max_log_a_entry);
	gtk_widget_show (max_log_cd8_entry);

	gtk_widget_show_all( main_window );

	/* since the button is present, start this off! */
	theState->writeOn = 0;

	/* these set here since widgets realized by now */
	gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (followTopStrains), theState->followTopStrains != 0);
	gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (plotTopStrains), theState->plotTopStrains != 0);
	gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (plotFirstStrains), theState->plotFirstStrains != 0);
	gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (writeEnable), theState->writeOn != 0);
	gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (scrollAxes), theState->scrollAxes != 0);
	gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (plot_vt_button), theState->plotVt != 0);
	gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (plot_act_button), theState->plotAct != 0);
	gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (plot_CD8_button), theState->plotCD8s != 0);
	gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (show_params_button), theState->show_params != 0);
	gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (show_stats_button), theState->show_stats != 0);
	gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (colorByHamming), theState->colorByHamming != 0);

	if (theState->followTopStrains == 0 && theState->plotTopStrains == 0 && theState->plotFirstStrains == 0)
	{
	    gtk_widget_hide( GTK_WIDGET(top_strain_box) );
	    gtk_widget_hide( GTK_WIDGET(first_strain_box) );
	}
	if (theState->plotFirstStrains == 0)
	    gtk_widget_hide( GTK_WIDGET(first_strain_box) );
}

/**************************************************************************
 *  * The following section contains utility function definitions.
 *   **************************************************************************/


#endif

void renderText(double x, double y, string instring)
{
  /*
   * Show font description string.
   */
#ifndef NO_GUI
  glRasterPos2f ((float)x, (float)y);
  glListBase (font_list_base);
  glCallLists (instring.length(), GL_UNSIGNED_BYTE, instring.c_str());
#else
  glRasterPos2f ((float)x, (float)y);
  float rpos[4];
  glGetFloatv(GL_CURRENT_RASTER_POSITION ,rpos);

  print (our_font, rpos[0], rpos[1], instring.c_str());
#endif
}

void renderGraphText(double x, double y, string instring)
{
  /*
   * Show font description string.
   */
#ifndef NO_GUI
  glRasterPos2f ((float)x, (float)y);
  glListBase (font_list_base);
  glCallLists (instring.length(), GL_UNSIGNED_BYTE, instring.c_str());
#else
  glRasterPos2f ((float)x, (float)y);
  float rpos[4];
  glGetFloatv(GL_CURRENT_RASTER_POSITION ,rpos);

  print (our_font, rpos[0]+40., rpos[1]+40., instring.c_str());
#endif
}

void draw_filled_box(GLfloat theColor[],double x, double y, double w, double h)
{

  glBegin(GL_POLYGON);
  glColor4fv(theColor);
  glVertex2d(x,y);
  glVertex2d(x,y+h);
  glVertex2d(x,y+h);
  glVertex2d(x+w,y+h);
  glVertex2d(x+w,y+h);
  glVertex2d(x+w,y);
  glVertex2d(x+w,y);
  glVertex2d(x,y);
  glEnd();
  // draw outline
  glLineWidth(2);
  glBegin(GL_LINES);
  glColor4fv(Black);
  glVertex2d(x,y);
  glVertex2d(x,y+h);
  glVertex2d(x,y+h);
  glVertex2d(x+w,y+h);
  glVertex2d(x+w,y+h);
  glVertex2d(x+w,y);
  glVertex2d(x+w,y);
  glVertex2d(x,y);
  glEnd();
}

void draw_outline()
{

  // draw outline
  glLineWidth(2);
  glBegin(GL_LINES);
  glColor4fv(Black);
  glVertex2d(-60.,-60.0);
  glVertex2d(-60.,60.0);
  glVertex2d(-60.,60.0);
  glVertex2d(60.,60.0);
  glVertex2d(60.,60.0);
  glVertex2d(60.,-60.0);
  glVertex2d(60.,-60.0);
  glVertex2d(-60.,-60.0);
  glEnd();
}

void draw_axes(int which)
{
  // displayed interval
  double start_time;

  double max_log_v_value;
  double min_log_v_value;

  double max_log_a_value;
  double min_log_a_value;

  double max_log_cd8_value;
  double min_log_cd8_value;

  double blood_factor=1;

  if (theState->displayCompartment == 0)
      blood_factor = theState->vol / 1000;

  // draw X axis line
  glLineWidth(2);
  glBegin(GL_LINES);
  glColor4fv(Black);
  glVertex2d(ORIGIN_X,ORIGIN_Y);
  glVertex2d(MAX_X_COORD,ORIGIN_Y);

  // draw Y axis line - 2
  glVertex2d(ORIGIN_X,ORIGIN_Y);
  glVertex2d(ORIGIN_X,MAX_Y_COORD);
  glVertex2d(MAX_X_COORD,ORIGIN_Y);
  glVertex2d(MAX_X_COORD,MAX_Y_COORD);
  glEnd();

  // draw 9 horiz division lines (1/10th of X max apart)
#ifdef DRAW_HORIZ_LINES
  for (int i=0; i < 10; i++) 
  {
      glLineWidth(2);
      glBegin(GL_LINES);
      glColor4fv(Black);
      glVertex2d(ORIGIN_X,ORIGIN_Y+(int)(i*(1.0/theState->y1_ticks)*(MAX_Y_COORD-ORIGIN_Y)));
      glVertex2d(MAX_X_COORD,ORIGIN_Y+(int)(i*(1.0/theState->y1_ticks)*(MAX_Y_COORD-ORIGIN_Y)));
      glEnd();
  }
#endif

  //glColor4fv(Black);
  //renderGraphText(10,-10, "Time (days)");

  // are we showing the whole run?
  if (theState->plot_span == 0.0)
  {
	start_time = 0.0;
  }
  // if not, which interval are we showing?
  else
  {
	if (theState->scrollAxes)
	{
	    //always show last n days where n=theState->plot_span
	    //
	    if (theState->time+theState->time_bias <= theState->plot_span)
	    {
		start_time = 0.0;
	    }
	    else
	    {
		start_time = (theState->time+theState->time_bias)-theState->plot_span;
		start_time=MAX(0.0, start_time+theState->plot_bias);
	    }

	}
	else
	{
	    //divide current time by interval to get start time
	    //
	    if (theState->time+theState->time_bias > theState->plot_span)
	    {
		start_time = theState->plot_span * (double)((int)((theState->time+theState->time_bias)/theState->plot_span));
		start_time=MAX(0.0, start_time+theState->plot_bias);
	    }
	    else
	    {
		start_time = 0.0;
	    }

	}
  }
  double tick_interval = (theState->plot_span)/((double)(theState->x_ticks));
  double first_tick = (double)(((int)(start_time / tick_interval))*tick_interval);
  if (start_time > first_tick)
	first_tick+=tick_interval;

  double startx = ORIGIN_X + ((first_tick-start_time)/theState->plot_span) * (MAX_X_COORD - ORIGIN_X);
  double delta_xcoord = (tick_interval/theState->plot_span) * (MAX_X_COORD - ORIGIN_X);
  double xcoord;
  for (int i=0; i < theState->x_ticks+1; i++) {
      char val_str[100];

      /* determine the next 'tick' based on start time */
      sprintf(val_str,"%g",first_tick+i*tick_interval);

      xcoord = startx + delta_xcoord*i;

      //xcoord = (100/(double)(theState->x_ticks+1))*(first_tick - start_time)/tick_interval+(100/(double)(theState->x_ticks+1))*i;
      //xcoord = ORIGIN_X + delta_xcoord*i;
      if (xcoord >= ORIGIN_X && xcoord <= MAX_X_COORD)
      {
	  glColor4fv(Black);
	  renderGraphText(xcoord, -5.0, val_str);
	  if (xcoord > ORIGIN_X && xcoord < MAX_X_COORD)
	  {
	      glLineWidth(2);
	      glBegin(GL_LINES);
	      glColor4fv(Black);
	      glVertex2d(xcoord,ORIGIN_Y);
	      glVertex2d(xcoord,ORIGIN_Y-2);
	      glEnd();
	  }
      }
  }
  if (theState->ART_start >= start_time && theState->ART_eff > 0 &&
	theState->ART_start <= start_time + theState->plot_span)
  {
      double xcoord = ((theState->ART_start-start_time)/theState->plot_span) * (MAX_X_COORD - ORIGIN_X);
      glColor4fv(DarkRed);
      renderGraphText(xcoord, -10.0, "ART ON");
      if (xcoord > ORIGIN_X && xcoord < MAX_X_COORD)
      {
	  glLineWidth(2);
	  glBegin(GL_LINES);
	  glColor4fv(DarkRed);
	  glVertex2d(xcoord,ORIGIN_Y);
	  glVertex2d(xcoord,MAX_Y_COORD);
	  glEnd();
      }
  }
  if (theState->ART_stop >= start_time && theState->ART_eff > 0 &&
	theState->ART_stop <= start_time + theState->plot_span)
  {
      double xcoord = ((theState->ART_stop-start_time)/theState->plot_span) * (MAX_X_COORD - ORIGIN_X);
      glColor4fv(Black);
      renderGraphText(xcoord, -10.0, "ART OFF");
      if (xcoord > ORIGIN_X && xcoord < MAX_X_COORD)
      {
	  glLineWidth(2);
	  glBegin(GL_LINES);
	  glColor4fv(Black);
	  glVertex2d(xcoord,ORIGIN_Y);
	  glVertex2d(xcoord,MAX_Y_COORD);
	  glEnd();
      }
  }
  if (theState->max_log_v_value > 0)
	max_log_v_value=theState->max_log_v_value;
  else
	max_log_v_value=(double)((int)(log10((double)theState->max_vl/blood_factor)+0.5));

  if (theState->min_log_v_value > 0)
	min_log_v_value=theState->min_log_v_value;
  else
	min_log_v_value=0;

  if (theState->max_log_a_value > 0)
	max_log_a_value=theState->max_log_a_value;
  else
	max_log_a_value=(double)((int)(log10((double)theState->max_act/blood_factor)+0.5));

  if (theState->min_log_a_value > 0)
	min_log_a_value=theState->min_log_a_value;
  else
	min_log_a_value=0;

  if (theState->max_log_cd8_value > 0)
	max_log_cd8_value=theState->max_log_cd8_value;
  else
	max_log_cd8_value=(double)((int)(log10((double)theState->max_cd8s/blood_factor)+0.5));

  if (theState->min_log_cd8_value > 0)
	min_log_cd8_value=theState->min_log_cd8_value;
  else
	min_log_cd8_value=0;

  // Take care of left axis labeling and title
  if (which == VT_GRAPH)
  {
      char title_str[100];
      glColor4fv(Black);
	  sprintf(title_str,"log 10 cell free HIV RNA copies by vol (1 ml)");

      renderGraphText(ORIGIN_X,90, title_str);
      renderGraphText(ORIGIN_X-12, MAX_Y_COORD-ORIGIN_Y+8, "Log");
      if (theState->followTopStrains || theState->plotTopStrains)
      {
	  int numTopStrains;
	  if (theState->displayCompartment == 0)
	      numTopStrains=theState->numTopStrains;
	  else
	      numTopStrains=theState->numTopStrains2;

	  char val_str[100];
	  if (theState->followTopStrains)
	      sprintf(val_str,"Following Strains ever in Top %d by VL (%d total of %d)",
		theState->maxTopStrains,numTopStrains,theState->numStrains);
	  else
	      sprintf(val_str,"Showing Current Top %d strains",theState->maxTopStrains);
	  renderGraphText(ORIGIN_X+5,85, val_str);
      }

      for (int i=0; i < theState->y1_ticks+1; i++)
      {
	  char val_str[100];
	  sprintf(val_str,"%g",min_log_v_value + ((max_log_v_value-min_log_v_value)/((double)(theState->y1_ticks)))*i);
	  renderGraphText(-8, (int)(i*(1.0/theState->y1_ticks)*(MAX_Y_COORD-ORIGIN_Y)), val_str);
      }
      for (int i=0; i < theState->y1_ticks+1; i++)
      {
	  glLineWidth(2);
	  glBegin(GL_LINES);
	  glColor4fv(Black);
	  glVertex2d(ORIGIN_X,ORIGIN_Y+(int)(i*(1.0/theState->y1_ticks)*(MAX_Y_COORD-ORIGIN_Y)));
	  glVertex2d(ORIGIN_X-2,ORIGIN_Y+(int)(i*(1.0/theState->y1_ticks)*(MAX_Y_COORD-ORIGIN_Y)));
	  glEnd();
      }
  }
  else if (which == ACT_GRAPH)
  {
      char title_str[100];
      glColor4fv(DarkBlue);
	  sprintf(title_str,"log 10 actively infected cells by vol (1 ml)");

      renderGraphText(ORIGIN_X,90, title_str);
      renderGraphText(ORIGIN_X-12, MAX_Y_COORD+8, "Log");

      if (theState->followTopStrains || theState->plotTopStrains)
      {
	  int numTopStrains;
	  if (theState->displayCompartment == 0)
	      numTopStrains=theState->numTopStrains;
	  else
	      numTopStrains=theState->numTopStrains2;

	  char val_str[100];
	  if (theState->followTopStrains)
	      sprintf(val_str,"Following Strains ever in Top %d by VL (%d total of %d)",
		theState->maxTopStrains,numTopStrains,theState->numStrains);
	  else
	      sprintf(val_str,"Showing Current Top %d strains",theState->maxTopStrains);
	  renderGraphText(ORIGIN_X+5,85, val_str);
      }

      for (int i=0; i < theState->y2_ticks+1; i++)
      {
	  char val_str[100];
	  sprintf(val_str,"%g",min_log_a_value + ((max_log_a_value-min_log_a_value)/((double)(theState->y2_ticks)))*i);
	  renderGraphText(-8, (int)(i*(1.0/theState->y2_ticks)*(MAX_Y_COORD-ORIGIN_Y)), val_str);
      }
      for (int i=0; i < theState->y2_ticks+1; i++)
      {
	  glLineWidth(2);
	  glBegin(GL_LINES);
	  glColor4fv(Black);
	  glVertex2d(ORIGIN_X,ORIGIN_Y+(int)(i*(1.0/theState->y2_ticks)*(MAX_Y_COORD-ORIGIN_Y)));
	  glVertex2d(ORIGIN_X-2,ORIGIN_Y+(int)(i*(1.0/theState->y2_ticks)*(MAX_Y_COORD-ORIGIN_Y)));
	  glEnd();
      }
  }
  else if (which == CD8_GRAPH)
  {
      char val_str[100];
      glColor4fv(Black);
      renderGraphText(ORIGIN_X,90, string("Log CD8 Cells"));
      for (int i=0; i < theState->y3_ticks+1; i++)
      {
	sprintf(val_str,"%g",min_log_cd8_value + ((max_log_cd8_value-min_log_cd8_value)/((double)(theState->y3_ticks)))*i);
	glColor4fv(Black);
	renderGraphText(ORIGIN_X-16, (int)(i*(1.0/theState->y3_ticks)*(MAX_Y_COORD-ORIGIN_Y)), val_str);
	glLineWidth(2);
	glBegin(GL_LINES);
	glColor4fv(Black);
	glVertex2d(ORIGIN_X,ORIGIN_Y+(int)(i*(1.0/theState->y3_ticks)*(MAX_Y_COORD-ORIGIN_Y)));
	glVertex2d(ORIGIN_X-2,ORIGIN_Y+(int)(i*(1.0/theState->y3_ticks)*(MAX_Y_COORD-ORIGIN_Y)));
	glEnd();
      }
      if(theState->immun_model == 3)
      {
	renderGraphText(MAX_X_COORD+8,90, string("Log"));
	if (theState->points != NULL)
	{
	    sprintf(val_str,"Total=%d",
		theState->points->cd8[MAX(theState->points->valid-1,0)]);
	    renderGraphText(ORIGIN_X,85,val_str);
	}
	for (int i=0; i < theState->y3_ticks+1; i++)
	{
	  sprintf(val_str,"%d",i);
	  glColor4fv(Black);
	  renderGraphText(MAX_X_COORD+4, (int)(i*(1.0/theState->y3_ticks)*(MAX_Y_COORD-ORIGIN_Y)), val_str);
	  glLineWidth(2);
	  glBegin(GL_LINES);
	  glColor4fv(Black);
	  glVertex2d(MAX_X_COORD,ORIGIN_Y+(int)(i*(1.0/theState->y3_ticks)*(MAX_Y_COORD-ORIGIN_Y)));
	  glVertex2d(MAX_X_COORD+2,ORIGIN_Y+(int)(i*(1.0/theState->y3_ticks)*(MAX_Y_COORD-ORIGIN_Y)));
	  glEnd();
	}
      }
  }

}
void show_params()
{
	char val_str[100];
	glColor4fv(DarkBlue);
	renderGraphText(ORIGIN_X,90, string("Simulation Parameters"));
	int ycoord=85;
	int ydelta=5;

	sprintf(val_str,"Beta mean: %e", theState->Bt_mean);
	renderGraphText(ORIGIN_X+10,ycoord, val_str);
	ycoord-=ydelta;
	if(theState->Bt_k > 0)
	{
	    sprintf(val_str,"Beta k: %lf", theState->Bt_k);
	    renderGraphText(ORIGIN_X+10,ycoord, val_str);
	    ycoord-=ydelta;
	}
	if(theState->Bt_max > 0)
	{
	    sprintf(val_str,"Beta max: %e", theState->Bt_max);
	    renderGraphText(ORIGIN_X+10,ycoord, val_str);
	    ycoord-=ydelta;
	}
	if(theState->lifespan_decay > 0)
	{
	    sprintf(val_str,"lifespan decay: %lf", theState->lifespan_decay);
	    renderGraphText(ORIGIN_X+10,ycoord, val_str);
	    ycoord-=ydelta;
	}
	sprintf(val_str,"Junk DNA fract: %lf", theState->junk);
	renderGraphText(ORIGIN_X+10,ycoord, val_str);
	ycoord-=ydelta;
	sprintf(val_str,"Mutation rate: %lf", theState->mu);
	renderGraphText(ORIGIN_X+10,ycoord, val_str);
	ycoord-=ydelta;
	sprintf(val_str,"ART efficiency: %e", theState->ART_eff);
	renderGraphText(ORIGIN_X+10,ycoord, val_str);
	ycoord-=ydelta;
	sprintf(val_str,"Production rate: %lf", theState->pi);
	renderGraphText(ORIGIN_X+10,ycoord, val_str);
	ycoord-=ydelta;
	sprintf(val_str,"Gamma: %lf", theState->gam);
	renderGraphText(ORIGIN_X+10,ycoord, val_str);
	ycoord-=ydelta;
	if(theState->eclipse > 0)
	{
	    sprintf(val_str,"eclipse: %lf", theState->eclipse);
	    renderGraphText(ORIGIN_X+10,ycoord, val_str);
	    ycoord-=ydelta;
	}
	if(theState->dA_scale > 0)
	{
	    sprintf(val_str,"dA scale: %lf", theState->dA_scale);
	    renderGraphText(ORIGIN_X+10,ycoord, val_str);
	    ycoord-=ydelta;
	}
	if(theState->dA_ic50_0 > 0)
	{
	    sprintf(val_str,"dA IC50: %lf", theState->dA_ic50_0);
	    renderGraphText(ORIGIN_X+10,ycoord, val_str);
	    ycoord-=ydelta;
	}

	if(theState->immun_model > 0)
	{
	    sprintf(val_str,"Theta: %lf", theState->theta);
	    renderGraphText(ORIGIN_X+10,ycoord, val_str);
	    ycoord-=ydelta;
	    sprintf(val_str,"Delta: %lf", theState->delta);
	    renderGraphText(ORIGIN_X+10,ycoord, val_str);
	    ycoord-=ydelta;
	    sprintf(val_str,"Immune IC50: %lf", theState->dImmun_IC50_0);
	    renderGraphText(ORIGIN_X+10,ycoord, val_str);
	    ycoord-=ydelta;

	}
	switch (theState->immun_model)
	{
	    case 0:
		if (theState->dA_ic50 > 0)
		    renderGraphText(ORIGIN_X+10,ycoord, "Strain-wise Immune Pressure via dA");
		else
		    renderGraphText(ORIGIN_X+10,ycoord, "No Immune Model");
		break;
	    case 1:
		renderGraphText(ORIGIN_X+10,ycoord, "Global Immune Model");
		break;
	    case 2:
		renderGraphText(ORIGIN_X+10,ycoord, "Strain-wise Immune Model");
		break;
	    case 3:
		renderGraphText(ORIGIN_X+10,ycoord, "Combined Immune Model");
		break;
	    case 4:
		renderGraphText(ORIGIN_X+10,ycoord, "Strain-group Immune Model by HD");
		break;
	    case 5:
		renderGraphText(ORIGIN_X+10,ycoord, "Strain-group Immune Model by birth time");
		break;
	    default:
		renderGraphText(ORIGIN_X+10,ycoord, "Unknown Immune Model");
		break;
	}
	ycoord-=ydelta;

	if(theState->immun_model == 1 || theState->immun_model == 3)
	{
	    sprintf(val_str,"Global Immune killing: %lf", theState->dImmun);
	    renderGraphText(ORIGIN_X+10,ycoord, val_str);
	    ycoord-=ydelta;
	}
	if(theState->immun_model == 2 || theState->immun_model >= 3)
	{
	    sprintf(val_str,"Strain-wise Immune killing: %lf", theState->dImmun_s);
	    renderGraphText(ORIGIN_X+10,ycoord, val_str);
	    ycoord-=ydelta;
	}
	if(theState->immun_model != 0)
	{
	    sprintf(val_str,"Intial CD8 Mean: %lf", theState->cd8_mean);
	    renderGraphText(ORIGIN_X+10,ycoord, val_str);
	    ycoord-=ydelta;
	    if(theState->cd8_k > 0)
	    {
		sprintf(val_str,"Intial CD8 k: %lf", theState->cd8_k);
		renderGraphText(ORIGIN_X+10,ycoord, val_str);
		ycoord-=ydelta;
	    }
	    if(theState->use_cd8s)
	    {
		sprintf(val_str,"CD8 wax ic50: %lf", theState->dImmun_IC50_0);
		renderGraphText(ORIGIN_X+10,ycoord, val_str);
		ycoord-=ydelta;

		sprintf(val_str,"CD8 Kappa: %lf", theState->kappa);
		renderGraphText(ORIGIN_X+10,ycoord, val_str);
		ycoord-=ydelta;

		sprintf(val_str,"CD8 wane IC50: %lf", theState->cd8_ic50_0);
		renderGraphText(ORIGIN_X+10,ycoord, val_str);
		ycoord-=ydelta;
	    }
	}
}
void show_stats()
{
	char val_str[100];
	glColor4fv(DarkBlue);
	renderGraphText(ORIGIN_X,90, string("Simulation Outputs"));
	int ycoord=85;
	int ydelta=5;

        if (theState->points == NULL || theState->points->valid == 0) 
	    return;

	sprintf(val_str,"Viable Strains: %d", theState->numStrains);
	renderGraphText(ORIGIN_X+10,ycoord, val_str);
	ycoord-=ydelta;

	double blood_factor=1;

	if (theState->displayCompartment == 0)
	    blood_factor = theState->vol / 1000;

	double val;
	val = (theState->points->vt_vi[MAX(theState->points->valid-1,0)] >0)?
		log10(theState->points->vt_vi[MAX(theState->points->valid-1,0)]/blood_factor):0;
	sprintf(val_str,"Log Viable VL: %lf", val);
	renderGraphText(ORIGIN_X+10,ycoord, val_str);
	ycoord-=ydelta;

	val = (theState->points->vt_non[MAX(theState->points->valid-1,0)] >0)?
		log10(theState->points->vt_non[MAX(theState->points->valid-1,0)]/blood_factor):0;
	sprintf(val_str,"Log Non-viable VL: %lf", val);
	renderGraphText(ORIGIN_X+10,ycoord, val_str);
	ycoord-=ydelta;

	val = (theState->points->junk_L0[MAX(theState->points->valid-1,0)] >0)?
		log10(theState->points->junk_L0[MAX(theState->points->valid-1,0)]/blood_factor):0;
	sprintf(val_str,"Log junk L0 DNA cells: %lf", val);
	renderGraphText(ORIGIN_X+10,ycoord, val_str);
	ycoord-=ydelta;

	val = (theState->points->junk_L1[MAX(theState->points->valid-1,0)] >0)?
		log10(theState->points->junk_L1[MAX(theState->points->valid-1,0)]/blood_factor):0;

	sprintf(val_str,"Log junk L1 DNA cells: %lf", val);
	renderGraphText(ORIGIN_X+10,ycoord, val_str);
	ycoord-=ydelta;

	val = (theState->points->junk_L2[MAX(theState->points->valid-1,0)] >0)?
		log10(theState->points->junk_L2[MAX(theState->points->valid-1,0)]/blood_factor):0;
	sprintf(val_str,"Log junk L2 DNA cells: %lf", val);
	renderGraphText(ORIGIN_X+10,ycoord, val_str);
	ycoord-=ydelta;

	val = (theState->points->act[MAX(theState->points->valid-1,0)] >0)?
		log10(theState->points->act[MAX(theState->points->valid-1,0)]/blood_factor):0;
	sprintf(val_str,"Log act cells: %lf", val);
	renderGraphText(ORIGIN_X+10,ycoord, val_str);
	ycoord-=ydelta;

	val = (theState->points->lat[MAX(theState->points->valid-1,0)] >0)?
		log10(theState->points->lat[MAX(theState->points->valid-1,0)]/blood_factor):0;
	sprintf(val_str,"Log latent cells: %lf", val);
	renderGraphText(ORIGIN_X+10,ycoord, val_str);
	ycoord-=ydelta;
}
//! [7]
void draw_graph(int which)
{
  double timeCoord, timeCoordp;

  double vtCoord;
  double vtNonCoord;
  double vtCoords[MAX_FOLLOW_STRAINS];
  double vtTarget;

  double cd8Coord;
  double cd8Coord2;
  double cd8Coords[MAX_FOLLOW_STRAINS];

  double actCoord;
  double actCoords[MAX_FOLLOW_STRAINS];

  // displayed interval (based on index of samples)
  double start_time;
  double end_time;
  double time_spread;

  int start_sample;
  int end_sample;
  int curr_sample;

  int k;

  double max_log_v_value;
  double min_log_v_value;
  double max_log_a_value;
  double min_log_a_value;
  double max_log_cd8_value;
  double min_log_cd8_value;
  double max_time;

  double blood_factor = theState->vol / 1000;

  if (which == SHOW_PARAMS)
  {
	show_params();
  }
  if (which == SHOW_STATS)
  {
	show_stats();
  }
  // No points yet? why bother...
  if (theState->points == NULL || theState->points->valid == 0) 
	return;

  if (theState->max_log_v_value > 0)
	max_log_v_value=theState->max_log_v_value;
  else
	max_log_v_value=(double)((int)(log10((double)theState->max_vl/blood_factor)+0.5));

  if (theState->min_log_v_value > 0)
	min_log_v_value=theState->min_log_v_value;
  else
	min_log_v_value=0;

  if (theState->min_log_a_value > 0)
	min_log_a_value=theState->min_log_a_value;
  else
	min_log_a_value=0;

  if (theState->max_log_a_value > 0)
	max_log_a_value=theState->max_log_a_value;
  else
	max_log_a_value=(double)((int)(log10((double)theState->max_act/blood_factor)+0.5));

  if (theState->min_log_cd8_value > 0)
	min_log_cd8_value=theState->min_log_cd8_value;
  else
	min_log_cd8_value=0;

  if (theState->max_log_cd8_value > 0)
	max_log_cd8_value=theState->max_log_cd8_value;
  else
	max_log_cd8_value=(double)((int)(log10((double)theState->max_cd8s/blood_factor)+0.5));


  start_sample=0;
  end_sample=MAX(theState->points->valid-1,0);
  curr_sample=MAX(theState->points->valid-1,0);
  time_spread = theState->max_time;

  // are we showing the whole run?
  // if not, which interval are we showing?
  if (theState->plot_span > 0.0)
  {
        time_spread = theState->plot_span;
	if (theState->scrollAxes)
	{
	    //always show last n days where n=theState->plot_span
	    //
	    if (theState->time+theState->time_bias <= theState->plot_span)
	    {
		start_time = 0.0;
		start_sample = 0;
	    }
	    else
	    {
		start_time = (theState->time+theState->time_bias)-theState->plot_span;
		start_time=MAX(0.0, start_time+theState->plot_bias);
		for (int j=0; j < theState->points->valid; j++) 
		{
		    if (theState->points->time[j] >= start_time)
		    {
			start_sample = j;
			break;
		    }
		}
	    }
	}
	else
	{
	    //divide current time by interval to get start time
	    //
	    if (theState->time+theState->time_bias > theState->plot_span)
	    {
		start_time = theState->plot_span * (double)((int)((theState->time+theState->time_bias)/theState->plot_span));
		start_time=MAX(0.0, start_time+theState->plot_bias);
		for (int j=0; j < theState->points->valid; j++) 
		{
		    if (theState->points->time[j] >= start_time)
		    {
			start_sample = j;
			break;
		    }
		}
	    }
	    else
	    {
		start_time = 0.0;
		start_sample = 0;
	    }
	}

	end_time = MIN(start_time + theState->plot_span,theState->time);

	for (int j=1; j < theState->points->valid; j++) 
	{
	    if (theState->points->time[j] >= theState->time+theState->time_bias)
	    {
		curr_sample = j-1;
		break;
	    }
	}

	for (int j=1; j < theState->points->valid; j++) 
	{
	    if (theState->points->time[j] >= end_time)
	    {
		end_sample = j-1;
		break;
	    }
	}
  }
  max_time = start_time + theState->plot_span;

    double val;
    glLineWidth(1);
    glBegin(GL_LINES);

    if (which == VT_GRAPH)
    {
      val = (theState->points->vt[start_sample] > 0)?log10((double)(theState->points->vt[start_sample]/blood_factor)):0;
      vtCoord = ((val-min_log_v_value)/(max_log_v_value-min_log_v_value)) *(MAX_Y_COORD-ORIGIN_Y);

      val = (theState->points->vt_non[start_sample] > 0)?log10((double)(theState->points->vt_non[start_sample]/blood_factor)):0;
      if (val < 0) val = 0;
      vtNonCoord = ((val-min_log_v_value)/(max_log_v_value-min_log_v_value)) *(MAX_Y_COORD-ORIGIN_Y);

      if (theState->followTopStrains || theState->plotTopStrains || theState->plotFirstStrains)
      {
	int numTopStrains;
	if (theState->displayCompartment == 0)
	    numTopStrains=theState->numTopStrains;
	else
	    numTopStrains=theState->numTopStrains2;

	for (k=0; k < numTopStrains;k++)
	{
	    if (theState->points->id_s[start_sample][k] >= 0)
	    {
		val = (theState->points->vt_s[start_sample][k] > 0)?log10((double)(theState->points->vt_s[start_sample][k]/blood_factor)):0;
	        if (val < 0) val = 0;
		vtCoords[k] = ((val-min_log_v_value)/(max_log_v_value-min_log_v_value)) *(MAX_Y_COORD-ORIGIN_Y);
	    }
	    else
		vtCoords[k] = -1;
	}
      }
    }
    else if (which == ACT_GRAPH)
    {
      val = (theState->points->act[start_sample] > 0)?log10((double)(theState->points->act[start_sample]/blood_factor)):0;
      if (val < 0) val = 0;
      actCoord = ((val-min_log_a_value)/(max_log_a_value-min_log_a_value)) *(MAX_Y_COORD-ORIGIN_Y);

      if (theState->followTopStrains || theState->plotTopStrains || theState->plotFirstStrains)
      {
	int numTopStrains;
	if (theState->displayCompartment == 0)
	    numTopStrains=theState->numTopStrains;
	else
	    numTopStrains=theState->numTopStrains2;
	for (k=0; k < numTopStrains;k++)
	{
	    if (theState->points->id_s[start_sample][k] >= 0)
		val = (theState->points->a_s[start_sample][k] > 0)?log10((double)(theState->points->a_s[start_sample][k]/blood_factor)):0;
	    else
		val = 0;
      	    if (val < 0) val = 0;
	    actCoords[k] = ((val-min_log_a_value)/(max_log_a_value-min_log_a_value)) *(MAX_Y_COORD-ORIGIN_Y);
	}
      }
    }
    else if (which == CD8_GRAPH)
    {
	if (theState->immun_model >= 1)
	{
	  val = (theState->points->cd8[start_sample] > 0)?log10((double)(theState->points->cd8[start_sample]/blood_factor)):0;
	  cd8Coord = ((val-min_log_cd8_value)/(max_log_cd8_value-min_log_cd8_value)) *(MAX_Y_COORD-ORIGIN_Y);
	}
	if(theState->immun_model >= 4)
	{
	  val = (theState->points->cd8_junk[start_sample] > 0)?log10((double)(theState->points->cd8_junk[start_sample]/blood_factor)):0;
	  cd8Coord2 = ((val-min_log_cd8_value)/(max_log_cd8_value-min_log_cd8_value)) *(MAX_Y_COORD-ORIGIN_Y);
	}
	if (theState->followTopStrains || theState->plotTopStrains || theState->plotFirstStrains || theState->immun_model >= 4)
	{
	    int numTopStrains;
	    if (theState->immun_model == 4)
		  numTopStrains=MIN(MAX_FOLLOW_STRAINS,theState->maxHammingDist/theState->cd8GroupSize);
	    else if (theState->immun_model == 5)
		  numTopStrains=MIN(theState->time/theState->cd8GroupSize,theState->tF/theState->cd8GroupSize);
	    if (theState->displayCompartment == 0)
		  numTopStrains=theState->numTopStrains;
	    else
		  numTopStrains=theState->numTopStrains2;
	    for (k=0; k < numTopStrains;k++)
	    {
		if (theState->points->id_s[start_sample][k] >= 0 || theState->immun_model >= 4)
		  val = (theState->points->cd8_s[start_sample][k] > 0)?log10((double)(theState->points->cd8_s[start_sample][k]/blood_factor)):0;
		else
		    val = 0;
		cd8Coords[k] = ((val-min_log_cd8_value)/(max_log_cd8_value-min_log_cd8_value)) *(MAX_Y_COORD-ORIGIN_Y);
	    }
	}
    }

    timeCoordp = ((double)(theState->points->time[start_sample]-start_time)/time_spread)
		    *MAX_X_COORD;
    for (int i = start_sample+1; i <= end_sample; i++) {
      timeCoord = ((double)(theState->points->time[i]-start_time)/time_spread) *MAX_X_COORD;

      if (which == VT_GRAPH)
      {
	  if (theState->plotTopStrains || theState->plotFirstStrains ||
		theState->followTopStrains)
	  {
	      int numTopStrains;
	      if (theState->displayCompartment == 0)
		  numTopStrains=theState->numTopStrains;
	      else
		  numTopStrains=theState->numTopStrains2;

	      for (k=0; k < numTopStrains;k++)
	      {
		  if (theState->colorByHamming && 
			theState->followTopStrains && theState->points->id_s[i][k] >=0){
		      double color_grad=theState->points->dh_s[i][k]/(double)theState->maxHammingDist;
		      int color_index = (int)(theState->max_strains * color_grad);
		      color_index=MIN(theState->max_strains-1,color_index);
		      glColor4fv(rainbow[color_index]);
		   }
		  else if (theState->followTopStrains && theState->points->id_s[i][k] >=0) {
		      glColor4fv(rainbow[theState->points->id_s[i][k]%theState->max_strains]);
		  }
		  else if (theState->points->color[i][k] >= 2 && theState->points->color[i][k]  < MAX_STRAIN_COLORS+2)
		      glColor4fv(pColors[theState->points->color[i][k]]);
		  else
		      glColor4fv(Black);
		  if (vtCoords[k] >=0 && theState->points->id_s[i][k] >= 0)
		  {
		      glVertex2d(timeCoordp,vtCoords[k]);
		      val = (theState->points->vt_s[i][k] > 0)?log10((double)(theState->points->vt_s[i][k]/blood_factor)):0;
		      if (val < 0) val = 0;
		      vtCoords[k] = ((val-min_log_v_value)/(max_log_v_value-min_log_v_value)) *(MAX_Y_COORD-ORIGIN_Y);
		      glVertex2d(timeCoord,vtCoords[k]);
		  }
		  else if (theState->points->id_s[i][k] >= 0)
		  {
		      val = (theState->points->vt_s[i][k] > 0)?log10((double)(theState->points->vt_s[i][k]/blood_factor)):0;
		      if (val < 0) val = 0;
		      vtCoords[k] = ((val-min_log_v_value)/(max_log_v_value-min_log_v_value)) *(MAX_Y_COORD-ORIGIN_Y);
		  }
	      }
	  }
	  {
	      /* write black line from last viral load reading to current one */
	      glColor4fv(Black);
	      glVertex2d(timeCoordp,vtCoord);
	      val = (theState->points->vt[i] > 0)?log10((double)(theState->points->vt[i]/blood_factor)):0;
	      if (val < 0) val = 0;
	      vtCoord = ((val-min_log_v_value)/(max_log_v_value-min_log_v_value)) *(MAX_Y_COORD-ORIGIN_Y);
	      glVertex2d(timeCoord,vtCoord);
	  }
   	  if (theState->showNonVia)
	  {
	      /* write gray line from last non-viable viral load reading to current one */
	      glColor4fv(Gray);
	      glVertex2d(timeCoordp,vtNonCoord);
	      val = (theState->points->vt_non[i] > 0)?log10((double)(theState->points->vt_non[i]/blood_factor)):0;
	      if (val < 0) val = 0;
	      vtNonCoord = ((val-min_log_v_value)/(max_log_v_value-min_log_v_value)) *(MAX_Y_COORD-ORIGIN_Y);
	      glVertex2d(timeCoord,vtNonCoord);
	  }
      }
      else if (which == ACT_GRAPH)
      {
	  if (theState->plotTopStrains || theState->plotFirstStrains ||
		theState->followTopStrains)
	  {
	      int numTopStrains;
	      if (theState->displayCompartment == 0)
		  numTopStrains=theState->numTopStrains;
	      else
		  numTopStrains=theState->numTopStrains2;

	      for (k=0; k < numTopStrains;k++)
	      {
		  if (theState->colorByHamming && 
			theState->followTopStrains && theState->points->id_s[i][k] >=0){
		      double color_grad=theState->points->dh_s[i][k]/(double)theState->maxHammingDist;
		      int color_index = (int)(theState->max_strains * color_grad);
		      color_index=MIN(theState->max_strains-1,color_index);
		      glColor4fv(rainbow[color_index]);
		   }
		  else if (theState->followTopStrains && theState->points->id_s[i][k] >=0)
		      glColor4fv(rainbow[theState->points->id_s[i][k]%theState->max_strains]);
		  else if (theState->points->color[i][k] >= 2 && theState->points->color[i][k]  < MAX_STRAIN_COLORS+2)
		      glColor4fv(pColors[theState->points->color[i][k]]);
		  else
		      glColor4fv(Black);
		  glVertex2d(timeCoordp,actCoords[k]);
		  if (theState->points->id_s[i][k] >= 0)
		      val = (theState->points->a_s[i][k] > 0)?log10((double)(theState->points->a_s[i][k]/blood_factor)):0;
		  else
		    val = 0;
		  if (val < 0) val = 0;
		  actCoords[k] = ((val-min_log_a_value)/(max_log_a_value-min_log_a_value)) *(MAX_Y_COORD-ORIGIN_Y);
		  glVertex2d(timeCoord,actCoords[k]);
	      }
	  }
	  {
	      /* write dark red line from last viral load reading to current one */
	      glColor4fv(Black);
	      glVertex2d(timeCoordp,actCoord);
	      val = (theState->points->act[i] > 0)?log10((double)(theState->points->act[i]/blood_factor)):0;
	      if (val < 0) val = 0;
	      actCoord = ((val-min_log_a_value)/(max_log_a_value-min_log_a_value)) *(MAX_Y_COORD-ORIGIN_Y);
	      glVertex2d(timeCoord,actCoord);
	  }
      }
      else if (which == CD8_GRAPH)
      {
	  if (theState->immun_model == 1)
	  {
	      /* write dark red line from last viral load reading to current one */
	      glColor4fv(Black);
	      glVertex2d(timeCoordp,cd8Coord);
	      val = (theState->points->cd8[i] > 0)?log10((double)(theState->points->cd8[i]/blood_factor)):0;
	      cd8Coord = ((val-min_log_cd8_value)/(max_log_cd8_value-min_log_cd8_value)) *(MAX_Y_COORD-ORIGIN_Y);
	      glVertex2d(timeCoord,cd8Coord);
	  }
	  else if (theState->immun_model > 1)
	  {
	      if (theState->plotTopStrains || theState->plotFirstStrains ||
		theState->followTopStrains)
	      {
		  int numTopStrains;
		  if (theState->immun_model == 4)
		      numTopStrains=MIN(MAX_FOLLOW_STRAINS,theState->maxHammingDist/theState->cd8GroupSize);
		  else if (theState->immun_model == 5)
		      numTopStrains=MIN(theState->time/theState->cd8GroupSize,theState->tF/theState->cd8GroupSize);
		  else if (theState->displayCompartment == 0)
		      numTopStrains=theState->numTopStrains;
		  else
		      numTopStrains=theState->numTopStrains2;

		  for (k=0; k < numTopStrains;k++)
		  {
		      if (theState->immun_model == 4)
		      {
			  double color_grad=(k*theState->cd8GroupSize)/(double)theState->maxHammingDist;
			  int color_index = (int)(theState->max_strains * color_grad);
			  color_index=MIN(theState->max_strains-1,color_index);
			  glColor4fv(rainbow[color_index]);
		      }
		      else if (theState->immun_model == 5)
		      {
			  double color_grad=(k*theState->cd8GroupSize)/(double)theState->tF;
			  int color_index = (int)(theState->max_strains * color_grad);
			  color_index=MIN(theState->max_strains-1,color_index);
			  glColor4fv(rainbow[color_index]);
		      }
		      else if (theState->colorByHamming && 
			    theState->followTopStrains && theState->points->id_s[i][k] >=0){
			  double color_grad=theState->points->dh_s[i][k]/(double)theState->maxHammingDist;
			  int color_index = (int)(theState->max_strains * color_grad);
			  color_index=MIN(theState->max_strains-1,color_index);
			  glColor4fv(rainbow[color_index]);
		       }
		      else if (theState->followTopStrains && theState->points->id_s[i][k] >=0) {
			  glColor4fv(rainbow[theState->points->id_s[i][k]%theState->max_strains]);
		      }
		      else if (theState->points->color[i][k] >= 2 && theState->points->color[i][k]  < MAX_STRAIN_COLORS+2)
			  glColor4fv(pColors[theState->points->color[i][k]]);
		      else
			  glColor4fv(Black);

		      if (theState->immun_model >= 4 || theState->points->id_s[i][k]>=0)
			  val = (theState->points->cd8_s[i][k] > 0)?log10((double)(theState->points->cd8_s[i][k]/blood_factor)):0;
		      else
			val = 0;

		      glVertex2d(timeCoordp,cd8Coords[k]);
		      cd8Coords[k] = ((val-min_log_cd8_value)/(max_log_cd8_value-min_log_cd8_value)) *(MAX_Y_COORD-ORIGIN_Y);
		      glVertex2d(timeCoord,cd8Coords[k]);
		  }
	      }
	      /* write black line from last total CD8 reading to current one */
	    if (theState->immun_model >= 3)
	    {
		  glColor4fv(Black);
		  glVertex2d(timeCoordp,cd8Coord);
		  val = (theState->points->cd8[i] > 0)?log10((double)(theState->points->cd8[i]/blood_factor)):0;
		  cd8Coord = ((val-min_log_cd8_value)/(max_log_cd8_value-min_log_cd8_value)) *(MAX_Y_COORD-ORIGIN_Y);
		  glVertex2d(timeCoord,cd8Coord);
	    }
	      /* write gray line from last total junk CD8 reading to current one */
	    if (theState->immun_model >= 4 && theState->showNonVia)
	    {
		  glColor4fv(Gray);
		  glVertex2d(timeCoordp,cd8Coord2);
		  val = (theState->points->cd8_junk[i] > 0)?log10((double)(theState->points->cd8_junk[i]/blood_factor)):0;
		  cd8Coord2 = ((val-min_log_cd8_value)/(max_log_cd8_value-min_log_cd8_value)) *(MAX_Y_COORD-ORIGIN_Y);
		  glVertex2d(timeCoord,cd8Coord2);
	    }
	  }
      }
      timeCoordp = timeCoord;
    }
    glEnd();

    if (which == VT_GRAPH)
    {
	// go back and place ALL target triangles where needed 
	double tri_time=start_time;
	while (tri_time < max_time)
	{
	  timeCoord = ((tri_time-start_time)/(time_spread)) *MAX_X_COORD;
	  tri_time += theState->sampleInterval;
	  if (theState->plotVt && theState->vt_points[theState->disp_patient] != NULL && (theState->Crit_type & 1))
	  {
	      for (int d=0; d < theState->num_vt_points[theState->disp_patient];d++)
	      {
		  if (fabs(theState->vt_points[theState->disp_patient][d].time - tri_time) 
				  < 0.5*theState->sampleInterval)
		  {
		    val = 0;
		    if (theState->vt_points[theState->disp_patient][d].vt > 0)
		    {
		      val=theState->vt_points[theState->disp_patient][d].vt;
		      vtTarget = (val/max_log_v_value) *(MAX_Y_COORD-ORIGIN_Y);
		      glBegin(GL_TRIANGLES);
		      glColor4fv(DarkRed);
		      glVertex2d(timeCoord,vtTarget);
		      glVertex2d(timeCoord-1,vtTarget-2);
		      glVertex2d(timeCoord+1,vtTarget-2);
		      glEnd();
		   }
		  }
	      }
	  }
	}
    }

    // draw the time triangle
    if (curr_sample <= end_sample && curr_sample >= start_sample)
    {
      timeCoord = ((double)(theState->points->time[curr_sample]-start_time)/time_spread)
			    *MAX_X_COORD;
      glBegin(GL_TRIANGLES);
      glColor4fv(Black);
      glVertex2d(timeCoord,2);
      glVertex2d(timeCoord-1,0);
      glVertex2d(timeCoord+1,0);
      glEnd();
    }

    // label the traces (if strain-specific only)
    if (!legend)
    {
	glColor4fv(Black);
	char time_str[100];
	sprintf(time_str,"Time = %0.2f (days)",theState->time+theState->time_bias);
	renderGraphText(40,-15, time_str);

	char val_str[100];
	if (theState->plotTopStrains || theState->plotFirstStrains)
	{
	    int numTopStrains;
	    if (theState->displayCompartment == 0)
		  numTopStrains=theState->numTopStrains;
	    else
		  numTopStrains=theState->numTopStrains2;

	    for (int i=0; i < numTopStrains+2; i++)
	    {
		if (i == 0)
		{
		    glColor4fv(Black);
		    strcpy(val_str,"total");
		}
		else if (i == 1 && theState->showNonVia)
		{
		    glColor4fv(Gray);
		    strcpy(val_str,"non-viable");
		}
		else if (i > 1 && theState->points->id_s[curr_sample][i-2] >= 0 &&
			theState->points->color[curr_sample][i-2] >= 1 && theState->points->color[curr_sample][i-2]  < MAX_STRAIN_COLORS+2)
		{
		    glColor4fv(pColors[theState->points->color[curr_sample][i-2]]);
		    sprintf(val_str,"strain %d",theState->points->id_s[curr_sample][i-2]);
		}
		else
		    continue;

		renderGraphText(MAX_X_COORD+4, 2+MAX_Y_COORD - (int)(i*(1.0/theState->maxTopStrains)*(MAX_Y_COORD-ORIGIN_Y)), val_str);
		glBegin(GL_TRIANGLES);
		glVertex2d(MAX_X_COORD+4,MAX_Y_COORD - (int)(i*(1.0/theState->maxTopStrains)*(MAX_Y_COORD-ORIGIN_Y)));
		glVertex2d(MAX_X_COORD+16,MAX_Y_COORD - (int)(i*(1.0/theState->maxTopStrains)*(MAX_Y_COORD-ORIGIN_Y)));
		glVertex2d(MAX_X_COORD+16,-2+MAX_Y_COORD - (int)(i*(1.0/theState->maxTopStrains)*(MAX_Y_COORD-ORIGIN_Y)));
		glVertex2d(MAX_X_COORD+16,-2+MAX_Y_COORD - (int)(i*(1.0/theState->maxTopStrains)*(MAX_Y_COORD-ORIGIN_Y)));
		glVertex2d(MAX_X_COORD+4,-2+MAX_Y_COORD - (int)(i*(1.0/theState->maxTopStrains)*(MAX_Y_COORD-ORIGIN_Y)));
		glVertex2d(MAX_X_COORD+4,MAX_Y_COORD - (int)(i*(1.0/theState->maxTopStrains)*(MAX_Y_COORD-ORIGIN_Y)));
		glEnd();
	    }
	}
	else if (which == VT_GRAPH || which == CD8_GRAPH)
	{
	    if (theState->followTopStrains)
	    {
		for (int i=0; i < 2+theState->colorBoxes; i++)
		{
		    strcpy(val_str,"");
		    if (i > 1)
		    {
			glColor4fv(rainbow[(theState->max_strains*(i-2)/theState->colorBoxes)]);
			glBegin(GL_TRIANGLES);
			glVertex2d(MAX_X_COORD+16,-2+MAX_Y_COORD - (int)((double)i*(MAX_Y_COORD-ORIGIN_Y)/(1+theState->colorBoxes)));
			glVertex2d(MAX_X_COORD+4,-2+MAX_Y_COORD - (int)((double)i*(MAX_Y_COORD-ORIGIN_Y)/(1+theState->colorBoxes)));
			glVertex2d(MAX_X_COORD+4,MAX_Y_COORD - (int)((double)i*(MAX_Y_COORD-ORIGIN_Y)/(1+theState->colorBoxes)));
			glEnd();
			glColor4fv(rainbow[(theState->max_strains*(i-1)/theState->colorBoxes)-1]);
			glBegin(GL_TRIANGLES);
			glVertex2d(MAX_X_COORD+4,MAX_Y_COORD - (int)((double)i*(MAX_Y_COORD-ORIGIN_Y)/(1+theState->colorBoxes)));
			glVertex2d(MAX_X_COORD+16,MAX_Y_COORD - (int)((double)i*(MAX_Y_COORD-ORIGIN_Y)/(1+theState->colorBoxes)));
			glVertex2d(MAX_X_COORD+16,-2+MAX_Y_COORD - (int)((double)i*(MAX_Y_COORD-ORIGIN_Y)/(1+theState->colorBoxes)));
			glEnd();
		        if (theState->colorByHamming)
			    sprintf(val_str,"HD <%d",theState->maxHammingDist*(i-1)/theState->colorBoxes);
			else
			    sprintf(val_str,"<%d",theState->max_strains*(i-1)/theState->colorBoxes);
		    }
		    else if (i > 0) 
		    { 
			if (theState->showNonVia)
			{
			    glColor4fv(Gray);
			    glBegin(GL_TRIANGLES);
			    glVertex2d(MAX_X_COORD+4,MAX_Y_COORD - (int)((double)i*(MAX_Y_COORD-ORIGIN_Y)/(1+theState->colorBoxes)));
			    glVertex2d(MAX_X_COORD+16,MAX_Y_COORD - (int)((double)i*(MAX_Y_COORD-ORIGIN_Y)/(1+theState->colorBoxes)));
			    glVertex2d(MAX_X_COORD+16,-2+MAX_Y_COORD - (int)((double)i*(MAX_Y_COORD-ORIGIN_Y)/(1+theState->colorBoxes)));
			    glVertex2d(MAX_X_COORD+16,-2+MAX_Y_COORD - (int)((double)i*(MAX_Y_COORD-ORIGIN_Y)/(1+theState->colorBoxes)));
			    glVertex2d(MAX_X_COORD+4,-2+MAX_Y_COORD - (int)((double)i*(MAX_Y_COORD-ORIGIN_Y)/(1+theState->colorBoxes)));
			    glVertex2d(MAX_X_COORD+4,MAX_Y_COORD - (int)((double)i*(MAX_Y_COORD-ORIGIN_Y)/(1+theState->colorBoxes)));
			    glEnd();
			    strcpy(val_str,"V non-viable");
			}
			else 
			    strcpy(val_str,"Strain #");
		    } 
		    else if (i==0)
		    {
			glColor4fv(Black);
			glBegin(GL_TRIANGLES);
			glVertex2d(MAX_X_COORD+4,MAX_Y_COORD - (int)((double)i*(MAX_Y_COORD-ORIGIN_Y)/(1+theState->colorBoxes)));
			glVertex2d(MAX_X_COORD+16,MAX_Y_COORD - (int)((double)i*(MAX_Y_COORD-ORIGIN_Y)/(1+theState->colorBoxes)));
			glVertex2d(MAX_X_COORD+16,-2+MAX_Y_COORD - (int)((double)i*(MAX_Y_COORD-ORIGIN_Y)/(1+theState->colorBoxes)));
			glVertex2d(MAX_X_COORD+16,-2+MAX_Y_COORD - (int)((double)i*(MAX_Y_COORD-ORIGIN_Y)/(1+theState->colorBoxes)));
			glVertex2d(MAX_X_COORD+4,-2+MAX_Y_COORD - (int)((double)i*(MAX_Y_COORD-ORIGIN_Y)/(1+theState->colorBoxes)));
			glVertex2d(MAX_X_COORD+4,MAX_Y_COORD - (int)((double)i*(MAX_Y_COORD-ORIGIN_Y)/(1+theState->colorBoxes)));
			glEnd();
			strcpy(val_str,"V all strains");
		    }

		    renderGraphText(MAX_X_COORD+4, 2+MAX_Y_COORD - (int)((double)i*(MAX_Y_COORD-ORIGIN_Y)/(1+theState->colorBoxes)), val_str);
		}
	    }
	    else
	    {
		glColor4fv(Black);
		strcpy(val_str,"V total");
		renderGraphText(MAX_X_COORD+4, MAX_Y_COORD+2, val_str);
		glBegin(GL_TRIANGLES);
		glVertex2d(MAX_X_COORD+4,MAX_Y_COORD);
		glVertex2d(MAX_X_COORD+16,MAX_Y_COORD);
		glVertex2d(MAX_X_COORD+16,-2+MAX_Y_COORD);
		glVertex2d(MAX_X_COORD+16,-2+MAX_Y_COORD);
		glVertex2d(MAX_X_COORD+4,-2+MAX_Y_COORD);
		glVertex2d(MAX_X_COORD+4,MAX_Y_COORD);
		glEnd();

		if (theState->showNonVia)
		{
		    glColor4fv(Gray);
		    strcpy(val_str,"V non-viable");
		    renderGraphText(MAX_X_COORD+4, MAX_Y_COORD+2-(MAX_Y_COORD-ORIGIN_Y)/(1+theState->colorBoxes), val_str);
		    glBegin(GL_TRIANGLES);
		    glVertex2d(MAX_X_COORD+4,MAX_Y_COORD-(MAX_Y_COORD-ORIGIN_Y)/(1+theState->colorBoxes));
		    glVertex2d(MAX_X_COORD+16,MAX_Y_COORD-(MAX_Y_COORD-ORIGIN_Y)/(1+theState->colorBoxes));
		    glVertex2d(MAX_X_COORD+16,-2+MAX_Y_COORD-(MAX_Y_COORD-ORIGIN_Y)/(1+theState->colorBoxes));
		    glVertex2d(MAX_X_COORD+16,-2+MAX_Y_COORD-(MAX_Y_COORD-ORIGIN_Y)/(1+theState->colorBoxes));
		    glVertex2d(MAX_X_COORD+4,-2+MAX_Y_COORD-(MAX_Y_COORD-ORIGIN_Y)/(1+theState->colorBoxes));
		    glVertex2d(MAX_X_COORD+4,MAX_Y_COORD-(MAX_Y_COORD-ORIGIN_Y)/(1+theState->colorBoxes));
		    glEnd();
		}
	    }
	}
	legend = true;
    }

    glFlush();
}

static bool
draw_routine ( GLfloat , GLfloat )
{
  /*** OpenGL BEGIN ***/
  glClear (GL_COLOR_BUFFER_BIT);

  /* Set the foreground colour. */
  glColor3f(0.0,0.0,0.0);

  /* count num active plots (1-6) */

  int plotCnt = 0;
  if (theState->plotVt)
      plotCnt++;
  if (theState->plotCD8s)
      plotCnt++;
  if (theState->plotAct)
      plotCnt++;
  if (theState->show_params)
      plotCnt++;
  if (theState->show_stats)
      plotCnt++;

  int rows, cols;
  switch(plotCnt)
  {
      case 1:
	  cols=1;
	  rows=1;
	  break;
      case 2:
	  cols=2;
	  rows=1;
	  break;
      case 3:
	  cols=3;
	  rows=1;
	  break;
      case 4:
	  cols=2;
	  rows=2;
	  break;
      default:
	  cols=3;
	  rows=2;
	  break;
  }
  double xfactor = 1./(double)cols;;
  double yfactor = 1./(double)rows;;

  double xshift_size = 120.;
  double yshift_size = 120.;

  double xshift = 0.;
  double yshift = 0.;

  glClear(GL_COLOR_BUFFER_BIT);
  glLoadIdentity();
  glPushMatrix();
  glScaled(xfactor, yfactor,0.0);

  if (plotCnt > 3)
    yshift = 60.;

  if (cols > 2)
  {
    xshift = -120.;
  }
  else if (cols > 1)
    xshift = -60.;

  glTranslated(xshift, yshift,0.0);
  plotCnt = 0;
  drawnTime = 0;
  legend = false;
  if (theState->plotVt)
  {
	//draw_outline();
	plotCnt++;
	glPushMatrix();
	glTranslated(-40.0, -40.0,0.0);
	glPushMatrix();
	draw_axes(VT_GRAPH);
	draw_graph(VT_GRAPH);
	glPopMatrix();
	glPopMatrix();
	if (plotCnt != cols)
	{
	  xshift= xshift_size;
	  yshift= 0.0;
	}
	else
	{
	  xshift=-((cols-1)*xshift_size);
	  yshift=-yshift_size;
	}
	glTranslated(xshift, yshift,0.0);
  }
  if (theState->show_params)
  {
	//draw_outline();
	plotCnt++;
	glPushMatrix();
	glTranslated(-40.0, -40.0,0.0);
	glPushMatrix();
	draw_graph(SHOW_PARAMS);
	glPopMatrix();
	glPopMatrix();
	if (plotCnt != cols)
	{
	  xshift= xshift_size;
	  yshift= 0.0;
	}
	else
	{
	  xshift=-((cols-1)*xshift_size);
	  yshift=-yshift_size;
	}
	glTranslated(xshift, yshift,0.0);
  }
  if (theState->show_stats)
  {
	//draw_outline();
	plotCnt++;
	glPushMatrix();
	glTranslated(-40.0, -40.0,0.0);
	glPushMatrix();
	draw_graph(SHOW_STATS);
	glPopMatrix();
	glPopMatrix();
	if (plotCnt != cols)
	{
	  xshift= xshift_size;
	  yshift= 0.0;
	}
	else
	{
	  xshift=-((cols-1)*xshift_size);
	  yshift=-yshift_size;
	}
	glTranslated(xshift, yshift,0.0);
  }
  if (theState->plotAct)
  {
	//draw_outline();
	plotCnt++;
	glPushMatrix();
	glTranslated(-40.0, -40.0,0.0);
	glPushMatrix();
	draw_axes(ACT_GRAPH);
	draw_graph(ACT_GRAPH);
	glPopMatrix();
	glPopMatrix();
	if (plotCnt != cols)
	{
	  xshift= xshift_size;
	  yshift= 0.0;
	}
	else
	{
	  xshift=-((cols-1)*xshift_size);
	  yshift=-yshift_size;
	}
	glTranslated(xshift, yshift,0.0);
  }
  if (theState->plotCD8s)
  {
	//draw_outline();
	plotCnt++;
	glPushMatrix();
	glTranslated(-40.0, -40.0,0.0);
	glPushMatrix();
	draw_axes(CD8_GRAPH);
	draw_graph(CD8_GRAPH);
	glPopMatrix();
	glPopMatrix();
	if (plotCnt != cols)
	{
	  xshift= xshift_size;
	  yshift= 0.0;
	}
	else
	{
	  xshift=-((cols-1)*xshift_size);
	  yshift=-yshift_size;
	}
	glTranslated(xshift, yshift,0.0);
  }
  glPopMatrix();

  
  glFlush ();

  return true;
}

int gui_main (int argc, char *argv[], settings *vars)
{
    pColors[0] = White;
    pColors[1] = DarkRed;
    pColors[2] = DarkGreen;
    pColors[3] = DarkBlue;
    pColors[4] = Orange;
    pColors[5] = Gold;
    pColors[6] = DarkPurple;
    pColors[7] = Pink;
    pColors[8] = Yellow;
    pColors[9] = Brick;
    pColors[10] = Gold;
    pColors[11] = Blue;

    theState = vars;

    for (int i=0; i < theState->max_strains; i++)
    {
	rainbow[i] = (GLfloat *)malloc(4*sizeof(GLfloat));
	assignRainbow(rainbow[i],i,theState->max_strains);
    }

    /* Initialize GTK. */
    gtk_init (&argc, &argv);

    /* Initialize GtkGLExt. */
    gtk_gl_init (&argc, &argv);

    make_gui();

    gtk_main ();

    return 0;

}//end of gui_main

void update_points(settings *vars, bool *snap_this_frame, int snapnum,
	unsigned int vi_Vt, unsigned int nonvi_Vt, generic_list<strain *> *topQs,
	unsigned int acts, unsigned int lats, unsigned int cd8s, 
	unsigned int junk_L0s, unsigned int junk_L1s, unsigned int junk_L2s)
{ 
    unsigned int Vt = vi_Vt+nonvi_Vt;

    // reset at start of the run
    if (vars->time == 0)
    {
	vars->points->valid=0;
    }

    if (vars->points->valid < vars->points->max_points-1) {
	vars->points->time[vars->points->valid] = vars->time;
	vars->points->vt[vars->points->valid] = vi_Vt+nonvi_Vt;
	vars->points->vt_vi[vars->points->valid] = vi_Vt;
	vars->points->vt_non[vars->points->valid] = nonvi_Vt;
	vars->points->act[vars->points->valid] = acts;
	vars->points->lat[vars->points->valid] = lats;
	vars->points->cd8[vars->points->valid] = cd8s;
	vars->points->junk_L0[vars->points->valid] = junk_L0s;
	vars->points->junk_L1[vars->points->valid] = junk_L1s;
	vars->points->junk_L2[vars->points->valid] = junk_L2s;
	while (vars->immun_model > 0 && cd8s > (unsigned int)vars->max_cd8s)
	    vars->max_cd8s = 2* vars->max_cd8s;

	if (vars->immun_model >= 4)
	{
	    vars->points->cd8_junk[vars->points->valid] = 
	        vars->cd8_groups[vars->displayCompartment][0].cd8_cells;

	    int max_index;
	    if (vars->immun_model == 4)
		max_index = MIN(MAX_FOLLOW_STRAINS,theState->maxHammingDist/theState->cd8GroupSize);
	    else if (vars->immun_model == 5)
		max_index = MIN((int)(theState->time/theState->cd8GroupSize),(int)(theState->tF/theState->cd8GroupSize));

	    for (int j=0; j <max_index; j++)
	    {
		int cd8_index = j+1;
		vars->points->cd8_s[vars->points->valid][j] = theState->cd8_groups[vars->displayCompartment][cd8_index].cd8_cells;
		if (vars->points->cd8_s[vars->points->valid][j] > vars->max_cd8s)
		    vars->max_cd8s = 2* vars->max_cd8s;
	    }
	    for (int j=max_index; j <MAX_FOLLOW_STRAINS; j++)
		vars->points->cd8_s[vars->points->valid][j] = 0;
	}
	for (int j=0; j <MAX_FOLLOW_STRAINS; j++)
	{
	    strain *the_strain = topQs->get_elem(j);
	    if (the_strain!=NULL)
	    {
		vars->points->id_s[vars->points->valid][j] = the_strain->get_strain_index();
		vars->points->vt_s[vars->points->valid][j] = (the_strain->get_via_virus(vars->displayCompartment)+the_strain->get_junk_virus(vars->displayCompartment));
		if (vars->immun_model == 2 || vars->immun_model == 3)
		{
		    vars->points->cd8_s[vars->points->valid][j] = the_strain->get_cd8_cells(vars->displayCompartment);
		    if (vars->points->cd8_s[vars->points->valid][j] > vars->max_cd8s)
			vars->max_cd8s = 2* vars->max_cd8s;

		}
		vars->points->a_s[vars->points->valid][j] = the_strain->get_act_cells(vars->displayCompartment);
		vars->points->dh_s[vars->points->valid][j] = the_strain->get_founder_dist();
	    }
	    else
	    {
		vars->points->id_s[vars->points->valid][j] = -1;
		vars->points->vt_s[vars->points->valid][j] = 0;
		if (vars->immun_model == 2 || vars->immun_model == 3)
		    vars->points->cd8_s[vars->points->valid][j] = 0;
		vars->points->a_s[vars->points->valid][j] = 0;
		vars->points->dh_s[vars->points->valid][j] = 0;
	    }
	    vars->points->color[vars->points->valid][j] = j+1;
	}
	vars->points->valid++;

	if (vars->time > vars->max_time)
	    vars->max_time = 2* vars->max_time;
	if (Vt > vars->max_vl)
	    vars->max_vl = 2* vars->max_vl;

	if (image != NULL)
	{
	    if (vars->time >= vars->NextRefresh) {
		updateGL();
		vars->NextRefresh+=vars->refresh;
	    }
	    if (vars->AutoSnapshot &&  
		((snapnum == 0) ||
		(vars->time>=vars->NextSnap))) { 
		*snap_this_frame=true;  // done after graph updates!
		updateGL();
		vars->NextSnap+=vars->SnapshotInterval;
	    }
	}
	else
	{
	    if (vars->AutoSnapshot &&  
		((snapnum == 0) ||
		(vars->time>=vars->NextSnap))) { 
		*snap_this_frame=true;  // done after graph updates!
		updateGL();
		vars->NextSnap+=vars->SnapshotInterval;
	    }
	}
    }
    /* check for need to realloc point array */
    vars->points->checkForRealloc();
}

void check_for_pause(settings *vars, bool *snap_this_frame)
{
      /* check for graph updates atleast once (and whenever paused) */
    do {
	if (gtk_events_pending ())
	    gtk_main_iteration();
	if (*snap_this_frame)
	{
	    char snapfile[100];
	    vars->snapshot++;
	    sprintf(snapfile,
		"snapshot_%d.png",vars->snapshot);
	    take_screenshot( NULL, snapfile);
	    *snap_this_frame=false;
	}
    } while (vars->pauseFlag == 1);
}
