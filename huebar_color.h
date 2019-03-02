#include <math.h>

/**
 * Total number of available HUE colors.
 */
const int NUMBER_OF_HUE_COLORS = 360;

/**
 * Representation of the RGB colors.
 */
typedef struct RGB {
	/**
	 * Red channel value.
	 */
	unsigned char r;

	/**
	 * Green channel value.
	 */
	unsigned char g;

	/**
	 * Blue channel value.
	 */
	unsigned char b;
} RGB;

RGB hue2rgb[NUMBER_OF_HUE_COLORS];

/**
 *
 * @param idx
 * @param r Red channel value.
 * @param g Green channel value.
 * @param b Blue channel value.
 * @param I
 */
void set_hue2rgb_channels(int idx, int r, int g, int b, double I) {
	if (idx < 0) {
		return;
	}

	if (idx >= NUMBER_OF_HUE_COLORS) {
		return;
	}

	I = 1;

	hue2rgb[idx].r = r * I;
	hue2rgb[idx].g = g * I;
	hue2rgb[idx].b = b * I;
}

/**
 * Initialization of the HUE to RGB converter.
 */
void init_hue2rgb() {
	/**/
	int s1 = 0;

	/**/
	int s2 = NUMBER_OF_HUE_COLORS / 3;

	/**/
	int s3 = 2 * NUMBER_OF_HUE_COLORS / 3;

	/**/
	int s4 = NUMBER_OF_HUE_COLORS;

	/**/
	int ss = 1 + (NUMBER_OF_HUE_COLORS / 3);

	for(int i = 0; i <= ss/2; i++) {
		double a = sin( i * M_PI/ss );

		int S = 255 * a;

		double j = 2.0 * (double)i / (double)ss;

		set_hue2rgb_channels(s1+i, 255, S, 0, j);
		set_hue2rgb_channels(s2+i, 0, 255, S, j);
		set_hue2rgb_channels(s3+i, S, 0, 255, j);
		set_hue2rgb_channels(s2-i, S, 255, 0, j);
		set_hue2rgb_channels(s3-i, 0, S, 255, j);
		set_hue2rgb_channels(s4-i, 255, 0, S, j);
	}
}

/**
 * Creation of HUE bar.
 *
 * @param bar_length Length of the bar.
 *
 * @return Newly allocated bar.
 */
RGB* create_color_huebar(int bar_length) {
	//TODO Allocation should be done with a clear idea how the memory will be released.
	RGB *bar = (RGB *)malloc(sizeof(RGB) * bar_length);

	/**/
	int sx = bar_length;

	init_hue2rgb();

	for(int x=0; x < sx; x++) {
		int i = (x * NUMBER_OF_HUE_COLORS) / sx;
		int r = hue2rgb[i].r;
		int g = hue2rgb[i].g;
		int b = hue2rgb[i].b;

		bar[x].r = r;
		bar[x].g = g;
		bar[x].b = b;
	}

	return bar;
}
