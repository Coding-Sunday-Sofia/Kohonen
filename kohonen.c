#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <stdarg.h>
#include <math.h>
#include <time.h>
#include "huebar_color.h"

/* START Kohonen Algorithm defines and global variables. */

/** SOM width. */
#define MAP_WIDTH 80

/** SOM height. */
#define MAP_HEIGHT 60

/** Scaling factor. */
#define NORMALIZATION_VALUE 10000

/** Power of two macros. */
#define pow2(x) ((x) * (x))

/** Big integer nubmer. */
#define BIG_NUM          999999999999999999

/** File line length limit. */
#define LINE_SIZE        300

/** In C boolean values are not predefined. */
#define FALSE            0

/** In C boolean values are not predefined. */
#define TRUE             1

/** Alias of unsigned integer. */
#define uint             unsigned int

/** Alias of unsigned long. */
#define uint64           unsigned long

/** KNN neuron definition. */
typedef struct Neuron {
	/** List of neurons. */
	unsigned int* components;
} Neuron;

/** Best matching unit coordinates structure. */
typedef struct BMU {
	/** Column. */
	unsigned int x_coord;

	/** Row. */
	unsigned int y_coord;
} BMU;

/** Coordinates of the centroid in the best macthing unit cluster. */
typedef struct CentroidBMU {
	/** Column. */
	unsigned int x_coord;

	/** Row. */
	unsigned int y_coord;

	/** Nubmer of units in the cluster. */
	int count;
} CentroidBMU;

/** Cartesian coordinates structure. */
typedef struct Coordinate {
	/** X axis value. */
	float x;

	/** Y axis value. */
	float y;
} Coordinate;

/** Training samples structure. */
typedef struct Sample {
	/** List of training samples. */
	unsigned int* components;
} Sample;

/** Map as matrix of neurons. */
Neuron** map;

/** List of samples. */
Sample* samples;

/** */
int total_components;

/** */
char** components_name;

/** */
uint* samples_max_components_values;

/** */
uint* samples_min_components_values;

/** */
int initial_radius = 80;

/** */
float round_radius;

/* END Kohonen definitions and global variables. */

/* START k-Means Algorithm defines and global variables. */

/** Euclidean distance macros definition. */
#define distance(i, j) (datax(j) - datax(i)) * (datax(j) - datax(i)) + (datay(j) - datay(i)) * (datay(j) - datay(i))

/** Definition of a boolean type. */
typedef int bool;

/** Total number of training samples. */
int total_samples;

/* END k-Means definitions and global variables. */

/* START Kohonen algorithm methods. */

/**
 * Pseudo randum number generator scaling function.
 *
 * @param min Lower bound of the generated numbers.
 * @param max Upper bound of the generated numbers.
 *
 * @return Uniform distributed pseudo-random number into the specified range.
 */
#define randr(min, max) ((max - min + 1) * (double)rand() / RAND_MAX + min)

/**
 * Separate string in tokens.
 *
 * @param a_str String to split.
 * @param a_delim Splitting delimiters.
 * @param total_items Total number of expected tokens.
 *
 * @return List of tokens.
 */
char** explode_string(char* a_str, const char a_delim, int* total_items) {
	char** result    = 0;
	size_t count     = 0;
	char* tmp        = a_str;
	char* last_comma = 0;
	char delim[2];
	delim[0] = a_delim;
	delim[1] = 0;

	/* Count how many elements will be extracted. */
	while (*tmp) {
		if (a_delim == *tmp) {
			count++;
			last_comma = tmp;
		}
		tmp++;
	}

	/* Add space for trailing token. */
	count += last_comma < (a_str + strlen(a_str) - 1);
	*total_items = count;

	/*
	 * Add space for terminating null string so caller knows where the list of
	 * returned strings ends.
	 */
	count++;

	result = malloc(sizeof(char*) * count);

	if (result) {
		size_t idx  = 0;
		char* token = strtok(a_str, delim);

		while (token) {
			*(result + idx++) = strdup(token);
			token = strtok(0, delim);
		}
		*(result + idx) = 0;
	}

	return result;
}

/**
 * Convert string to integer.
 *
 * @param a String value.
 *
 * @return Integer value.
 */
int stringToInteger(char a[]) {
	int c, sign, offset, n;

	/* Handle negative integers. */
	if (a[0] == '-') {
		sign = -1;
	}

	/* Set starting position to convert. */
	if (sign == -1) {
		offset = 1;
	} else {
		offset = 0;
	}

	n = 0;

	for (c = offset; (a[c] != '\0') && (a[c] != '\n'); c++) {
		n = n * 10 + a[c] - '0';
	}

	if (sign == -1) {
		n = -n;
	}

	return n;
}

/**
 * Load training samples from a file.
 *
 * @param filename File name for the file with samples.
 */
void load_and_initialize_samples(char *filename) {
	FILE *file;
	char line[LINE_SIZE];
	char** line_components;
	int ch, x, y, z, total_line_components;
	uint value;

	total_samples = 0;

	file = fopen(filename, "r");
	while(!feof(file)) {
		ch = fgetc(file);
		if(ch == '\n') {
			total_samples++;
		}
	}
	total_samples--;
	fclose(file);

	printf("\n\nTotal samples: %d\n\n", total_samples);

	/*
	 * Variable 'samples' is the data structure used to store the sample points
	 * for Kohonen algorithm process.
	 */
	samples = (Sample *) malloc(sizeof(Sample) * total_samples);

	for(int i = 0; i < total_samples; i++) {
		samples[i].components = (unsigned int *) malloc(sizeof(unsigned int) * total_components);
	}

	/* Read data from file into array. */
	file = fopen(filename, "rt");
	fgets(line, LINE_SIZE, file);
	components_name = explode_string(line, ',', &total_components);

	samples_max_components_values = (uint*) malloc(sizeof(uint) * total_components);
	samples_min_components_values = (uint*) malloc(sizeof(uint) * total_components);

	/* Initialize max values array to blank. */
	for(int e = 0; e < total_components; e++) {
		samples_max_components_values[e] = 0;
		samples_min_components_values[e] = 0;
	}

	/* Load values from file. */
	for (int i = 0; i < total_samples; i++) {
		fgets(line, LINE_SIZE, file);
		line_components = explode_string(line, ',', &total_line_components);
		for(int e = 0; e < total_line_components; e++) {
			value = (unsigned int)(stringToInteger(line_components[e]));
			samples[i].components[e] = value;
			if(samples_max_components_values[e] < value) {
				samples_max_components_values[e] = value;
			}
			if(samples_min_components_values[e] > value) {
				samples_min_components_values[e] = value;
			}
		}
	}

	fclose(file);

	/* Normalize values based on max value of each component from 0 to 255. */
	for (int i = 0; i < total_samples; i++) {
		for(int e = 0; e < total_components; e++) {
			value = samples[i].components[e];
			samples[i].components[e] = (unsigned int)(((value - samples_min_components_values[e]) * NORMALIZATION_VALUE)/(samples_max_components_values[e] - samples_min_components_values[e]));
		}
	}
}

/**
 * SOM initialization.
 */
void initialize_som_map() {
	map = (Neuron **) malloc(sizeof(Neuron *) * MAP_WIDTH);

	int x, y;
	for(x = 0; x < MAP_WIDTH; x++) {
		map[x] = (Neuron *) malloc(sizeof(Neuron) * MAP_HEIGHT);
	}

	for(x = 0; x < MAP_WIDTH; x++) {
		for(y = 0; y < MAP_HEIGHT; y++) {
			map[x][y].components = (unsigned int *) malloc(sizeof(unsigned int) * total_components);

			for(int i=0; i<total_components; i++) {
				map[x][y].components[i] = randr(0,NORMALIZATION_VALUE);
			}
		}
	}
}

/**
 * Random sample selection.
 *
 * @return Pointer ot selected sample.
 */
Sample* pick_random_sample() {
	int i = randr(0, total_samples-1);
	return &samples[i];
}

/**
 * Index sample selection.
 *
 * @return Pointer ot selected sample.
 */
Sample* pick_sample(int i) {
	return &samples[i];
}

/**
 * Calculates distance between sample and SOM neuron.
 *
 * @param sample Sample pointer.
 * @param neuron Neuron pointer.
 *
 * @return Calculated distance.
 */
uint distance_between_sample_and_neuron(Sample *sample, Neuron *neuron) {
	unsigned int euclidean_distance = 0;
	unsigned int component_diff;

	for(int i = 0; i < total_components; i++) {
		component_diff = sample->components[i] - neuron->components[i];
		euclidean_distance += pow2(component_diff);
	}

	return euclidean_distance;
	//return sqrt(euclidean_distance);
}

/**
 * Searchs for best mathcing unit.
 *
 * @param sample Sample pointer.
 *
 * @return Pointer to the best matching unit.
 */
BMU* search_bmu(Sample *sample) {
	uint max_dist=999999999;
	uint dist = 0;
	BMU *bmu = (BMU *) malloc(sizeof(BMU));

	for(int x = 0; x < MAP_WIDTH; x++) {
		for(int y = 0; y < MAP_HEIGHT; y++) {
			dist = distance_between_sample_and_neuron(sample, &map[x][y]);
			if(dist < max_dist) {
				bmu->x_coord = x;
				bmu->y_coord = y;
				max_dist = dist;
			}
		}
	}

	return bmu;
}

/**
 * Calculates disttance between particular coordinates.
 *
 * @param p1 First point.
 * @param p2 Second point.
 *
 * @return Calculated distance.
 */
float get_coordinate_distance(Coordinate *p1, Coordinate *p2) {
	float x_sub = (p1->x) - (p2->x);
	float y_sub = (p1->y) - (p2->y);

	return sqrt(x_sub*x_sub + y_sub*y_sub);
}

/**
 * Creates coordinates object from given x and y values.
 *
 * @param x X value.
 * @param y Y value.
 *
 * @return Coordinates object.
 */
Coordinate* new_coordinate(float x, float y) {
	Coordinate *coordinate = malloc(sizeof(Coordinate));

	coordinate->x = x;
	coordinate->y = y;

	return coordinate;
}

/**
 * Select neuron at particualr positon.
 *
 * @param x X coordinate.
 * @param y Y coordinate.
 * @param sample Sample pointer.
 * @param scale Scaling factor.
 */
void scale_neuron_at_position(int x, int y, Sample *sample, double scale) {
	float neuron_prescaled, neuron_scaled;
	Neuron *neuron = &map[x][y];

	for(int i=0; i<total_components; i++) {
		neuron_prescaled = neuron->components[i] * (1.0f-scale);
		neuron_scaled = (sample->components[i] * scale) + neuron_prescaled;
		neuron->components[i] = (int)neuron_scaled;
	}
}

/**
 * Select neighbours.
 *
 * @param bmu Best matching unit pointer.
 * @param sample Sample pointer.
 * @param t
 */
void scale_neighbors(BMU *bmu, Sample *sample, float t) {
	float iteration_radius = roundf((float)(round_radius)*(1.0f-t));
	Coordinate *outer = new_coordinate(iteration_radius,iteration_radius);
	Coordinate *center = new_coordinate(0.0f,0.0f);
	float distance_normalized = get_coordinate_distance(center,outer);
	float distance;
	double scale;
	int x_coord;
	int y_coord;

	for(float y = -iteration_radius; y<iteration_radius; y++) {
		for(float x = -iteration_radius; x<iteration_radius; x++) {
			/* Out of bounds. */
			if(y + bmu->y_coord < 0) {
				continue;
			}

			/* Out of bounds. */
			if(y + bmu->y_coord >= MAP_HEIGHT) {
				continue;
			}

			/* Out of bounds. */
			if(x + bmu->x_coord < 0) {
				continue;
			}

			/* Out of bounds. */
			if(x + bmu->x_coord >= MAP_WIDTH) {
				continue;
			}

			outer->x = x;
			outer->y = y;

			distance = get_coordinate_distance(outer,center) / distance_normalized;

			/* Gaussian function. */
			scale = exp(-1.0f * (pow(distance, 2.0f)) / 0.15f);
			/* Exponential regulated cosine function. */
			//scale = cos(distance) / exp( fabs(distance) );
			/* Fading cosine function. */
			//scale = (-M_PI/2 <= distance && distance <= +M_PI/2) ? ( cos(distance) ) : ( cos(distance)/fabs(distance) );

			/* It is needed +1 is to avoid divide by 0's. */
			scale /= (t*4.0f + 1.0f);

			x_coord = bmu->x_coord + x;
			y_coord = bmu->y_coord + y;
			scale_neuron_at_position(x_coord, y_coord, sample, scale);
		}
	}

	free(outer);
	free(center);
}

/**
 * Release of the allocated memory.
 */
void free_allocated_memory() {
	for(int e = 0; e < total_components; e++) {
		free(components_name[e]);
	}
	free(components_name);

	for(int i = 0; i < total_samples; i++) {
		free(samples[i].components);
	}
	free(samples);

	for(int x = 0; x < MAP_WIDTH; x++) {
		for(int y = 0; y < MAP_HEIGHT; y++) {
			free(map[x][y].components);
		}
		free(map[x]);
	}
	free(map);

	free(samples_max_components_values);
	free(samples_min_components_values);
}

/**
 * Strings concatenation.
 *
 * @param s1 First string.
 * @param s2 Second string.
 *
 * @return Concatenated string.
 */
char* concat(const char *s1, const char *s2) {
	/* Plus one for the zero-terminator. */
	char *result = malloc(strlen(s1) + strlen(s2) + 1);

	strcpy(result, s1);
	strcat(result, s2);

	return result;
}

/**
 * Output generated as HTML file.
 *
 * @param final_bmus List of final best matching units.
 * @param auto_reload Flag for auto reload JavaScript code.
 */
void output_html(BMU *final_bmus, bool auto_reload) {
	bool found_bmu = FALSE;
	int diff_val;
	int max_val;
	int min_val;
	int x_val;
	int y_val;
	int z_val;
	int i;
	int value;
	int red;
	int green;
	int blue;
	uint x;
	uint y;
	float e;
	RGB *huebar = create_color_huebar(255);

	FILE *f = fopen("generated_kohonen_map.html", "w");
	if (f == NULL) {
		printf("Error opening file!\n");
		exit(1);
	}

	if(auto_reload) {
		fprintf(f, "<html><head><script>setTimeout(function(){ window.location.reload(1); }, 2500);</script></head><body>");
	} else {
		fprintf(f, "<html><head></head><body>");
	}

	fprintf(f, "<br/><h2>Components</h2>");

	for(int c = 0; c < total_components; c++) {
		fprintf(f, "<div style=\"width:500px;height:250px\">");
		fprintf(f, "<h3>%s</h3>", components_name[c]);
		fprintf(f, "<div style='position: absolute;'><table style='border-collapse: collapse;'>");
		/* Search the min and max values of each vector component. */
		max_val = 0;
		min_val = 9999999;
		for(y = 0; y < MAP_HEIGHT; y++) {
			for(x = 0; x < MAP_WIDTH; x++) {
				x_val = (int)((map[x][y].components[c] * (samples_max_components_values[c] - samples_min_components_values[c]))/NORMALIZATION_VALUE);

				if(min_val > x_val) {
					min_val = x_val;
				}

				if(max_val < x_val) {
					max_val = x_val;
				}
			}
		}

		diff_val = max_val - min_val;

		for(y = 0; y < MAP_HEIGHT; y++) {
			fprintf(f, "<tr>");
			for(x = 0; x < MAP_WIDTH; x++) {
				value = (int)((map[x][y].components[c] * (samples_max_components_values[c] - samples_min_components_values[c]))/NORMALIZATION_VALUE);
				x_val = (int)(((value - min_val) * 255)/diff_val);
				RGB *color = &huebar[x_val];
				fprintf(f, "<td style='width:3px;height:3px;background-color:rgb(%d,%d,%d);' title='%d'></td>", color->r, color->g, color->b, value);
			}
			fprintf(f, "</tr>");
		}
		fprintf(f, "</table></div>");

		fprintf(f, "<div style='position: absolute; left: 420px;'><table style='border-collapse: collapse;'>");
		for(e = 0.0f; e < 255.0f; e+=2.86f) {
			RGB *color = &huebar[(int)e];
			if(e == 0.0f) {
				fprintf(f, "<tr><td style='width:3px;height:1px;background-color:rgb(%d,%d,%d);'><div style=\"position:absolute; width:50px;top:0px;text-align:left;\">&nbsp; %d</div><div style=\"position: absolute; text-align: left; width: 50px; top: 85px;\">&nbsp; %d</div><div style=\"position:absolute; width:50px;bottom:0px;text-align:left;\">&nbsp; %d</div></td></tr>", color->r, color->g, color->b, min_val, min_val + (max_val-min_val)/2, max_val);
			} else {
				fprintf(f, "<tr><td style='width:3px;height:1px;background-color:rgb(%d,%d,%d);'></td></tr>", color->r, color->g, color->b);
			}
		}
		fprintf(f, "</table></div>");

		fprintf(f, "</div>");
	}

	fprintf(f, "</tr></table>");

	free(huebar);
	fprintf(f, "</body></html>");
	fclose(f);
}

/* END Kohonen algorithm methods. */

/**
 * Appcliation single entry point.
 *
 * @param argc Nubmer of command line arguments.
 * @param argv Command line arguments as list of strings.
 */
int main(int argc, char **argv) {
	char *filename = argv[1];
	load_and_initialize_samples(filename);

	int MAX_TRAINING_ROUNDS = 10;
	float ROUND_INC = 1.0f/(float)(MAX_TRAINING_ROUNDS);
	float r = 0.0f;

	/* Number of iterations is 30 times the number of input samples. */
	int MAX_ITER_PER_ROUND = total_samples * 30;
	float T_INC = 1.0f/(float)(MAX_ITER_PER_ROUND);
	float t = 0.0f;
	BMU *bmu;
	BMU *final_bmus;
	int iteration_num;
	int round_num;
	Sample *sample;

	/* Seed of pseudo-random number generator. */
	srand(time(NULL));

	final_bmus = (BMU *)malloc(sizeof(BMU) * total_samples);

	initialize_som_map();
	//output_html(final_bmus, TRUE);

	round_radius = initial_radius;
	round_num = 0;

	while(r < 1.0f) {
		round_radius = (r == 0.0f) ? initial_radius : (round_radius/2.0f);
		printf("\nROUND %d/%d | INITIAL RADIUS: %d | RADIUS: %f\n", round_num, MAX_TRAINING_ROUNDS, initial_radius, round_radius);

		iteration_num = 0;
		t = 0.0f;

		while(t < 1.0f) {
			sample = pick_random_sample();
			/* Best matching unit. */
			bmu = search_bmu(sample);
			scale_neighbors(bmu, sample, t);
			free(bmu);

			t += T_INC;
			iteration_num++;

			//output_html(final_bmus, TRUE);
			//usleep(100000);
		}

		r += ROUND_INC;
		round_num++;
	}

	/* Save the BMU coordinates in the SOM map. */
	for (int i = 0; i < total_samples; i++) {
		sample = pick_sample(i);
		/* Best matching unit. */
		bmu = search_bmu(sample);

		final_bmus[i].x_coord = (uint)bmu->x_coord;
		final_bmus[i].y_coord = (uint)bmu->y_coord;

		free(bmu);
	}

	output_html(final_bmus, FALSE);

	free(final_bmus);
	free_allocated_memory();

	return 0;
}
