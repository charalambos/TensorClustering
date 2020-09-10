////////////////////////////////////////////////////////////////////////////////////
// Copyright © Charalambos "Charis" Poullis, charalambos@poullis.org    	  //
// This work can only be used under an exclusive license of the author.           //
////////////////////////////////////////////////////////////////////////////////////
#ifndef __IMAGE_CPP__
#define __IMAGE_CPP__

#include "Image.h"
#include "Utilities.h"
#include "GeometryProcessing.h"
#include "ImageProcessing.h"
#include <exiv2/exiv2.hpp>

float Image::gamma_table[65536];
float Image::current_gamma = -1.0f;
unsigned short Image::lookup_table[65536];
bool Image::already_computed = false;

Image::Image(bool _verbose)	{

	verbose = _verbose;
	data = 0x00;
	width = height = 0;
	type = FT_UNSUPPORTED_FORMAT;
	file_name = "";

	if (!already_computed) {
		///Initialize the gamma table
		initializeGammaTable(1.0f);
		///Initialize the lookup table
		initializeLookupTable();

		already_computed = !already_computed;
	}
}

Image::Image(Image const &img, bool _verbose)	{

	verbose = _verbose;

	width = img.getWidth();
	height = img.getHeight();

	///Allocate and initialize memory
	data = new Color*[height];
	for (int y=0;y<height;y++)	{
		data[y] = new Color[width];
		for (int x=0;x<width;x++)	{
			data[y][x] = img.getPixel(x,y);
		}
	}

	type = img.getFileType();
	file_name = "";

	if (!already_computed) {
		///Initialize the gamma table
		initializeGammaTable(1.0f);
		///Initialize the lookup table
		initializeLookupTable();

		already_computed = !already_computed;
	}

}

Image::Image(int _width, int _height, bool _verbose)	{

	verbose = _verbose;

	///Set the dimensions
	width = _width;
	height = _height;

	///Allocate and initialize memory
	data = new Color*[height];
	for (int y=0;y<height;y++)	{
		data[y] = new Color[width];
		for (int x=0;x<width;x++)	{
			data[y][x] = Color(0.0f,0.0f,0.0f,1.0f);
		}
	}

	type = FT_UNSUPPORTED_FORMAT;
	file_name = "";

	if (!already_computed) {
		///Initialize the gamma table
		initializeGammaTable(1.0f);
		///Initialize the lookup table
		initializeLookupTable();

		already_computed = !already_computed;
	}
}

Image::Image(int _width, int _height, float r, float g, float b, float a, bool _verbose)	{

	verbose = _verbose;
	///Set the dimensions
	width = _width;
	height = _height;

	///Allocate and initialize memory
	data = new Color*[height];
	for (int y=0;y<height;y++)	{
		data[y] = new Color[width];
		for (int x=0;x<width;x++)	{
			data[y][x] = Color(r,g,b,a);
		}
	}

	type = FT_UNSUPPORTED_FORMAT;
	file_name = "";

	if (!already_computed) {
		///Initialize the gamma table
		initializeGammaTable(1.0f);
		///Initialize the lookup table
		initializeLookupTable();

		already_computed = !already_computed;
	}
}

Image::~Image()	{
	///Erase anything previously existing
	erase();
}

bool Image::loadImage(std::string const &_file_name)	{

	file_name = _file_name;
	///Separate the file name first
	std::string name;
	FILE_TYPE type;
	if (!separate(file_name, name, type))	{
		return false;
	}

	if (data)	{
		///Erase anything previously existing
		erase();
	}

	if (type != FT_PFM) {
        ///Load the contents
        try {
            Magick::Image img;
            img.read(file_name);
            width = img.columns();
            height = img.rows();

            ///Copy over to the internal structure
            Magick::Quantum *pixel_cache = img.getPixels(0, 0, width, height);
            size_t channels = img.channels();
            data = new Color *[height];
            for (int y = 0; y < height; y++) {
                data[y] = new Color[width];
                for (int x = 0; x < width; x++) {
                    Magick::ColorRGB color = img.pixelColor(x, height - 1 - y);

                    float r = color.red();
                    float g = color.green();
                    float b = color.blue();
                    float a = 1.0f;
                    if (channels == 4) a = color.alpha();

                    data[y][x] = Color(r, g, b, a);
                }
            }

        }
        catch (Magick::Exception &e) {
            std::cout << "Loading image exception: " << e.what() << std::endl;
            return false;
        }
    } else{
        if (!loadPFM(file_name)) return false;
	}

	return true;
}

bool Image::saveImage(std::string const &_file_name)	{

	file_name = _file_name;
	///Separate the file name first
	std::string name;
	FILE_TYPE type;
	if (!separate(file_name, name, type))	{
		return false;
	}

	if (type != FT_PFM)	{
		try {
			Magick::Image img( Magick::Geometry(size_t(width), size_t(height)), Magick::ColorRGB(0.0,0.0,0.0,0.0));

			///Copy over to the imagemagick's structure
			for (int y=0;y<height;y++)	{
				for (int x=0;x<width;x++)	{
					float r = data[y][x].r();
					float g = data[y][x].g();
					float b = data[y][x].b();
					float a = data[y][x].a();

					img.pixelColor(x, height-1-y,
								   Magick::ColorRGB(r,g,b,a));

				}
			}
			img.write(file_name);
			type = FT_PNG;
		}
		catch( Magick::Exception &e)	{
			std::cout << "Saving image  exception: " << e.what() << std::endl;
			return false;
		}

	}
	else	{
		if (!savePFM(file_name)) return false;
	}

	return true;
}


int Image::getWidth() const	{

	return width;
}

int Image::getHeight() const	{

	return height;
}

std::string Image::getFileName()	{

	return file_name;
}

FILE_TYPE Image::getFileType() const	{

	return type;
}

void Image::setFileType(FILE_TYPE _type)	{

	type = _type;
	return;
}

bool Image::loadPFM(std::string const &file_name)	{

	bool result = false;

	FILE *infile = 0x00;
	infile = fopen(file_name.c_str(), "rb");
	if (!infile)	return false;

	unsigned int a = fgetc(infile);
	int b = fgetc(infile);
	int junk = fgetc(infile);
	float *val = new float[3];
	float *cval = new float[1];

	if ((a != 'P') || ((b != 'F') && (b != 'f')))
		goto LoadPFM_error_end;

	b = (b == 'F');		// 'F' = RGB,  'f' = monochrome

	int x, y;
	fscanf(infile, "%d %d%c", &x, &y, &junk);
	if ((x <= 0) || (y <= 0))
		goto LoadPFM_error_end;

	float scalefactor;
	fscanf(infile, "%f%c", &scalefactor, &junk);

	width = x;
	height = y;
	data = new Color*[height];
	for (int i=0;i<height;i++)	{
		data[i] = new Color[width];
	}

	if (!data)	{
		result = false;
		goto LoadPFM_error_end;
	}

	for (y = 0; y < height; y++)	{
		for (x = 0; x < width; x++)	{
			data[y][x].a() = 1.0f;
			if (b)			// color
			{
				if (fread(val, sizeof(float), 3, infile) != 3)	goto LoadPFM_error_end;
				if (scalefactor > 0.0f)		// MSB
				{
					REVERSEBYTES(&(val[0]), &data[y][x].r());
					REVERSEBYTES(&(val[1]), &data[y][x].g());
					REVERSEBYTES(&(val[2]), &data[y][x].b());
				}
				else						// LSB
				{
					data[y][x].r() = val[0];
					data[y][x].g() = val[1];
					data[y][x].b() = val[2];
				}
// // 				delete [] val;
			}
			else			// black and white
			{
				if (fread(cval, sizeof(float), 1, infile) != 1) goto LoadPFM_error_end;
				if (scalefactor > 0.0f)		// MSB
				{
					REVERSEBYTES(&(cval[0]), &data[y][x].r());
				}
				else						// LSB
					data[y][x].r() = cval[0];

				data[y][x].g() = data[y][x].r();
				data[y][x].b() = data[y][x].r();
// 				delete [] cval;
			}

		}
	}
	result = true;

LoadPFM_error_end:
	fclose(infile);

	delete [] val;
	delete [] cval;

	type = FT_PFM;

	return result;
}

bool Image::savePFM(std::string const &file_name)	{
	FILE *outfile = fopen(file_name.c_str(), "wb");
	if (!outfile)	return 0;

	fputc('P', outfile);
	fputc('F', outfile);
	fputc(0x0a, outfile);

	fprintf(outfile, "%d %d", width, height);
	fputc(0x0a, outfile);

	fprintf(outfile, "%f", -1.0f);
	fputc(0x0a, outfile);

	bool result = true;
	float *val = new float[3];
	for (int y = 0; y < height; y++)	{
		for (int x=0;x<width; x++)	{
			val[0] = data[y][x].r();
			val[1] = data[y][x].g();
			val[2] = data[y][x].b();
			fwrite(val, sizeof(float), 3, outfile);
		}
	}
	delete [] val;

	fclose(outfile);

	return true;
}

Color Image::getPixel(int x, int y) const	{
	if (outOfBounds(const_cast<Image *>(this),x,y))	{
		if (verbose) {
			std::cout << _format("Indices (%d,%d) out of bounds. Image size (%d,%d)", x, y, this->width, this->height)
					  << std::endl;
		}
		return Color(0.0f,0.0f,0.0f,1.0f);
	}
	return data[y][x];
}

///The range of values is between 0-1
void Image::setPixel(int x, int y, Color const &color)	{
	if (outOfBounds(const_cast<Image *>(this),x,y))	{
			if (verbose) {
				std::cout
						<< _format("Indices (%d,%d) out of bounds. Image size (%d,%d)", x, y, this->width, this->height)
						<< std::endl;
			}
			return;
	}
	data[y][x] = color;
	return;
}

void Image::perPixelNormalization(float _range)	{
    for (int y=0;y<height;y++)	{
		for (int x=0;x<width;x++)	{
			Vector3f val = color2vector3(getPixel(x,y));
			val.normalize();

			setPixel(x,y,Color(val(0),val(1),val(2),1.0f));
		}
	}

	return;
}

void Image::normalize(float _range)	{

	Color *min_max_val = range();

	///Compute the max and min values for each channel
	Color min_val = min_max_val[0];
	Color max_val = min_max_val[1];

	///Normalize each pixel between 0-1
	Color max_minus_min = max_val - min_val;
	for (int y=0;y<height;y++)	{
		for (int x=0;x<width;x++)	{
			Color val = getPixel(x,y);
			Color norm_val = _range * (val - min_val)/max_minus_min;
			///Reset the alpha to 1
			norm_val.a() = 1.0f;
			setPixel(x,y,norm_val);
		}
	}

	delete [] min_max_val;

	return;
}

void Image::normalize(Color *_range)	{
	Color *min_max_val = _range;
	for (int y=0;y<height;y++)	{
		for (int x=0;x<width;x++)	{
			Color val = getPixel(x,y);
			val.r() = min_max_val[0].r() + val.r()*(min_max_val[1].r()-min_max_val[0].r());
			val.g() = min_max_val[0].g() + val.g()*(min_max_val[1].g()-min_max_val[0].g());
			val.b() = min_max_val[0].b() + val.b()*(min_max_val[1].b()-min_max_val[0].b());
			setPixel(x,y,val);
		}
	}

	return;
}

void Image::clear(Color const &val)	{
	for (int y=0;y<height;y++)	{
		for (int x=0;x<width;x++)	{
			setPixel(x,y,val);
		}
	}
	return;
}

void Image::clear(float _r, float _g, float _b)	{
	clear(Color(_r, _g, _b));
	return;
}

void Image::add(Image *_img)	{
	if (width != _img->getWidth() || height != _img->getHeight()) return;

	for (int y=0;y<height;y++)	{
		for (int x=0;x<width;x++)	{
			Color val = getPixel(x,y) + _img->getPixel(x,y);
			setPixel(x,y,val);
		}
	}

	return;
}

void Image::add(int x, int y, Color const &val)	{
	Color pix_val = getPixel(x,y);
	setPixel(x,y,pix_val+val);
	return;
}

void Image::copy(Image *_img)	{
	if (width != _img->getWidth() || height != _img->getHeight())	{
		///Erase anything previously existing
		erase();

		width = _img->getWidth();
		height = _img->getHeight();

		///Allocate and copy memory
		data = new Color*[height];
		for (int y=0;y<height;y++)	{
			data[y] = new Color[width];
			for (int x=0;x<width;x++)	{
				setPixel(x,y,_img->getPixel(x,y));
			}
		}
	}
	else	{
		///Just copy the data
		for (int y=0;y<height;y++)	{
			for (int x=0;x<width;x++)	{
				setPixel(x,y,_img->getPixel(x,y));
			}
		}
	}

	return;
}

void normalizeSet(std::vector<Image *> &images, float range, bool ignore_value, Color value) {
        ///Normalizes a set of images based on the min and max of the entire set
        ///Find the max and min first
        Color max_value = Color(FLT_MIN,FLT_MIN,FLT_MIN);
        Color min_value = Color(FLT_MAX,FLT_MAX,FLT_MAX);
        for (int i=0;i<images.size();i++)      {
                for (int y=0;y<images[i]->getHeight();y++) {
                        for (int x=0;x<images[i]->getWidth();x++) {
                                if (ignore_value && images[i]->getPixel(x,y) == value) continue;
                                if (max_value.r() < images[i]->getPixel(x,y).r())  max_value.r() = images[i]->getPixel(x,y).r();
                                if (max_value.g() < images[i]->getPixel(x,y).g())  max_value.g() = images[i]->getPixel(x,y).g();
                                if (max_value.b() < images[i]->getPixel(x,y).b())  max_value.b() = images[i]->getPixel(x,y).b();
                                if (min_value.r() > images[i]->getPixel(x,y).r())  min_value.r() = images[i]->getPixel(x,y).r();
                                if (min_value.g() > images[i]->getPixel(x,y).g())  min_value.g() = images[i]->getPixel(x,y).g();
                                if (min_value.b() > images[i]->getPixel(x,y).b())  min_value.b() = images[i]->getPixel(x,y).b();
                        }
                }
        }

//      printf("max_value: %f %f %f\n",max_value.r,max_value.g,max_value.b);
//      printf("min_value: %f %f %f\n",min_value.r,min_value.g,min_value.b);
        //then normalize
        Color max_minus_min = max_value - min_value;
        if (max_minus_min.r() == 0.0f)    max_minus_min.r() = 1.0f;
        if (max_minus_min.g() == 0.0f)    max_minus_min.g() = 1.0f;
        if (max_minus_min.b() == 0.0f)    max_minus_min.b() = 1.0f;
        for (int i=0;i<images.size();i++)      {
                for (int y=0;y<images[i]->getHeight();y++) {
                        for (int x=0;x<images[i]->getWidth();x++) {
                                images[i]->setPixel(x, y, range * (images[i]->getPixel(x,y) - min_value)/max_minus_min);
                        }
                }
        }
        return;
}

void Image::initializeGammaTable(float gamma)	{
	if (current_gamma == gamma)	return;

	for (int x = 0; x < 65536; x++)	{
		gamma_table[x] = (float) pow((float) x / 65535.0f, gamma);
	}

	current_gamma = gamma;

	return;
}


void Image::initializeLookupTable()	{

	for (int x = 0; x < 65536; x++)	{
		float tval = (float)(::pow(short2float(x) , 1.0f / current_gamma) * 65535.0f);
		lookup_table[x] = (tval > 65535.0f ? (unsigned short) 65535.0f : (unsigned short) tval);
	}

	return;
}

Color Image::getPixelBilinear(float sx, float sy)	{

	int bx = (int) sx;	int by = (int) sy;		// integer part
	int bxo = bx + 1;	int byo = by + 1;
	sx -= bx; sy -= by;							// fractional part [0,1)

	bx %= width;								// slow, hopefully they are already positive
	by %= height;
	bxo %= width;
	byo %= height;

	float osx = 1.0f - sx;
	float osy = 1.0f - sy;

	return (data[by][bx]   * osx * osy +
		data[by][bxo]  * sx  * osy +
		data[byo][bx]  * osx * sy  +
		data[byo][bxo] * sx  * sy);
}

Color Image::getPixelBicubic(float sx, float sy)	{

	float A = -0.75f;
	Color val(0.0f,0.0f,0.0f,1.0f);
	float totalweight = 0.0f;

	int bx = (int) sx;
	int by = (int) sy;

	bx %= width;
	by %= height;

	for (int y = by - 1; y < by + 3; y++)	{
		float disty = fabs(sy - float(y));
		float yweight;
/*		if (disty <= 1)			// photoshop bicubic
			yweight = (( A + 2.0 ) * disty - ( A + 3.0 )) * disty * disty + 1.0;
		else
			yweight = (( A * disty - 5.0 * A ) * disty + 8.0 * A ) * disty - 4.0 * A;
*/
		if (disty <= 1)			// 'true' bicubic
			yweight  = ( ( disty - 9.0f/5.0f ) * disty -   1.0f/5.0f     ) * disty + 1.0f;
		else
			yweight  = ( ( -1.0f/3.0f * (disty-1) + 4.0f/5.0f     ) * (disty-1) -   7.0f/15.0f ) * (disty-1);

		for (int x = bx - 1; x < bx + 3; x++)	{
			float distx = fabs(sx - x);
			float weight;
/*			if (distx <= 1)		// photoshop bicubic
				xweight = (( A + 2.0 ) * distx - (A + 3.0)) * distx * distx +1.0;
			else
				xweight = (( A * distx - 5.0 * A ) * distx + 8.0 * A ) * distx - 4.0 * A;
*/
			if (distx <= 1.0f)		// 'true' bicubic
				weight  = (( ( distx - 9.0f/5.0f ) * distx -   1.0f/5.0f     ) * distx + 1.0f) * yweight;
			else
				weight  = (( ( -1.0f/3.0f * (distx-1) + 4.0f/5.0f     ) * (distx-1.0f) -   7.0f/15.0f ) * (distx-1.0f)) * yweight;

			val += data[y][x] * weight;
			totalweight += weight;
		}
	}
	val *= 1.0f / totalweight;
	return val;
}

bool outOfBounds(Image *img, int x, int y, int offset_x, int offset_y, int error_code) {

        int width = img->getWidth() - offset_x;
        int height = img->getHeight() - offset_y;
        if (x >= offset_x && x < width-offset_x &&
            y >= offset_y && y < height-offset_y) return false;
        if (error_code != -1) {
        		std::cout << _format("Indices (%d,%d) are out of bounds. Image size (%d,%d)", x, y, width, height)
						  << std::endl;
		}
        return true;
}

void Image::halfSize()	{
	///Create a temporary image
	Image *temp = new Image(width/2, height/2);

	for (int y = 0; y < temp->getHeight(); y++)	{
		int cy = 2 * y;
		for (int x = 0; x < temp->getWidth(); x++)	{
			int cx = 2 * x;
			temp->setPixel(x,y, (getPixel(cx,cy) + getPixel(cx,cy+1) + getPixel(cx+1,cy) + getPixel(cx+1,cy+1))*0.25f);
		}
	}

	///Copy the data over to this
	this->copy(temp);

	///Clean up
	delete temp;

	return;
}

void Image::doubleSize()	{
	Image *temp = new Image(width*2, height*2);

	for (int y = 0; y < temp->getHeight(); y++)	{
		float cy = (float(y) - 0.5f ) / 2.0f;
		for (int x = 0; x < temp->getWidth(); x++)	{
			float cx = (float(x) - 0.5f ) / 2.0f;
			temp->setPixel(x,y, getPixel(_round(cx), _round(cy)));
		}
	}


	///Copy the data over to this
	this->copy(temp);

	///Clean up
	delete temp;

	return;
}

void Image::doubleSizeBilinear()	{
	Image *temp = new Image(width*2, height*2);

	for (int y = 0; y < temp->getHeight(); y++)	{
		float cy = (float(y) - 0.5f ) / 2.0f;
		for (int x = 0; x < temp->getWidth(); x++)	{
			float cx = (float(x) - 0.5f ) / 2.0f;
			temp->setPixel(x,y, getPixelBilinear(cx, cy));
		}
	}


	///Copy the data over to this
	this->copy(temp);

	///Clean up
	delete temp;

	return;
}

void Image::binary(Color const &threshold)	{
	for (int y=0;y<height;y++)	{
		for (int x=0;x<width;x++)	{
			if (data[y][x].r() > threshold.r())	{
				data[y][x].r() = 1.0f;
			}
			else	{
				data[y][x].r() = 0.0f;
			}
			if (data[y][x].g() > threshold.g())	{
				data[y][x].g() = 1.0f;
			}
			else	{
				data[y][x].g() = 0.0f;
			}
			if (data[y][x].b() > threshold.b())	{
				data[y][x].b() = 1.0f;
			}
			else	{
				data[y][x].b() = 0.0f;
			}
		}
	}
	return;
}

void Image::scale(float f)	{
	for (int y=0;y<height;y++)	{
		for (int x=0;x<width;x++)	{
			data[y][x].r() *= f;
			data[y][x].g() *= f;
			data[y][x].b() *= f;
		}
	}
	return;
}

float Image::getPixelGrayscale(int x, int y, int index)	{
	if (index == -1 || index > 3)	{
		Color val = getPixel(x,y);
		return 0.2989f * val.r() + 0.5870f * val.g() + 0.1140f * val.b();
	}
	else	{
		Color val = getPixel(x,y);
		return val(index);
	}
}

void Image::grayscale(int index)	{
	for (int y=0;y<height;y++)	{
		for (int x=0;x<width;x++)	{
			float val = getPixelGrayscale(x,y, index);
			data[y][x] = Color(val);
		}
	}
	return;
}

Color *Image::range()	{
    ///Compute the max and min values for each channel
	Color max_val(-FLT_MAX,-FLT_MAX,-FLT_MAX);
	Color min_val(FLT_MAX,FLT_MAX,FLT_MAX);

	///Go through the image and find the max and min values
	for (int y=0;y<height;y++)	{
		for (int x=0;x<width;x++)	{
			Color val = getPixel(x,y);

			if (val(0) > max_val(0))	max_val(0) = val(0);
			if (val(1) > max_val(1))	max_val(1) = val(1);
			if (val(2) > max_val(2))	max_val(2) = val(2);

			if (val(0) < min_val(0))	min_val(0) = val(0);
			if (val(1) < min_val(1))	min_val(1) = val(1);
			if (val(2) < min_val(2))	min_val(2) = val(2);
		}
	}

	Color *min_max_val = new Color[2];
	min_max_val[0] = min_val;
	min_max_val[1] = max_val;
	return min_max_val;
}

void Image::blur(float variance)	{
	///Create the gaussian filter
	float constant = 1.0f/(2.0f*M_PI*variance);
	int size = 6 * int(sqrtf(variance));
	if (size%2==0)	size++;
	size = std::max(size, 3);
	Image *gaussian_kernel = new Image(size, size, 0.0f, 0.0f, 0.0f, 1.0f);
	for (int y=-size/2; y<=size/2; y++)	{
		for (int x=-size/2; x<=size/2; x++)	{
			float val = constant * exp( - (float(x*x) + float(y*y))/(2.0f*variance));
			gaussian_kernel->setPixel(x+size/2, y+size/2, Color(val, val, val));
		}
	}

	///Apply the filter
	Image *filtered_image = ImageProcessing::convolve(this, gaussian_kernel);


	this->copy(filtered_image);

	///Delete the kernel
	delete gaussian_kernel;

	return;
}

Image *Image::mathOp(Image *img, IMAGE_MATH_OP const &type)	{
	///Check if the size is the same
	if (this->getWidth() != img->getWidth() || this->getHeight() != img->getHeight())	return 0x00;

	Image *result = new Image(this->getWidth(), this->getHeight());
	for (int y=0;y<img->getHeight();y++)	{
		for (int x=0;x<img->getWidth();x++)	{
			Color val;
			switch(type)	{
			case IMAGE_MATH_OP_ADD:	val = this->getPixel(x,y) + img->getPixel(x,y);
									break;
			case IMAGE_MATH_OP_SUBTRACT:	val = this->getPixel(x,y) - img->getPixel(x,y);
											break;
			case IMAGE_MATH_OP_MULTIPLY:	val = this->getPixel(x,y) * img->getPixel(x,y);
											break;
			case IMAGE_MATH_OP_DIVIDE:		if (img->getPixel(x,y).r() != 0.0f)	{
												val.r() = this->getPixel(x,y).r() / img->getPixel(x,y).r();
											}
											if (img->getPixel(x,y).g() != 0.0f)	{
												val.g() = this->getPixel(x,y).g() / img->getPixel(x,y).g();
											}
											if (img->getPixel(x,y).b() != 0.0f)	{
												val.b() = this->getPixel(x,y).b() / img->getPixel(x,y).b();
											}
											break;
			}
			result->setPixel(x,y,val);
		}
	}

	return result;
}

Image *Image::negate()	{
	Image *result = new Image(this->getWidth(), this->getHeight());
	for (int y=0;y<this->getHeight();y++)	{
		for (int x=0;x<this->getWidth();x++)	{
			Color val = Color(0.0f,0.0f,0.0f)-this->getPixel(x,y);
			result->setPixel(x,y,val);
		}
	}
	return result;
}

void Image::RGBtoXYZ()	{
	for (int y=0;y<height;y++)	{
		for (int x=0;x<width;x++)	{
			float r = getPixel(x,y).r();
			float g = getPixel(x,y).g();
			float b = getPixel(x,y).b();

			if (r > 0.04045f){ r = pow((r + 0.055f) / 1.055f, 2.4f); }
			else { r = r / 12.92f; }
			if ( g > 0.04045f){ g = pow((g + 0.055f) / 1.055f, 2.4f); }
			else { g = g / 12.92f; }
			if (b > 0.04045){ b = pow((b + 0.055f) / 1.055f, 2.4f); }
			else {	b = b / 12.92; }

			r = r * 100.0f;
			g = g * 100.0f;
			b = b * 100.0f;

			//Observer. = 2°, Illuminant = D65

			float X = r * 0.4124f + g * 0.3576f + b * 0.1805f;
			float Y = r * 0.2126f + g * 0.7152f + b * 0.0722f;
			float Z = r * 0.0193f + g * 0.1192f + b * 0.9505f;

			setPixel(x,y,Color(X,Y,Z));
		}
	}

	return;
}

void Image::RGBtoLUV()	{
    ///convert first to XYZ
	RGBtoXYZ();

	///CIE chromaticity coordinates:
	float xn = 0.312713f;
	float yn = 0.329016f;
	///CIE luminance:
	float Yn = 1.0f;
	for (int y=0;y<height;y++)	{
		for (int x=0;x<width;x++)	{
			float X = getPixel(x,y).r();
			float Y = getPixel(x,y).g();
			float Z = getPixel(x,y).b();

			float un = 4.0f*xn / (-2.0f*xn + 12.0f*yn + 3.0f);
			float vn = 9.0f*yn / (-2.0f*xn + 12.0f*yn + 3.0f);
			float denom = (X + 15.0f*Y + 3.0f*Z);
			if (denom < TC_EPSILON)    denom = 1.0f;
			float u = 4.0f*X / denom;
			float v = 9.0f*Y / denom;
			float L = 116.0f * pow((Y/Yn),(1.0f/3.0f)) - 16.0f;
			float U = 13.0f*L*(u-un);
			float V = 13.0f*L*(v-vn);

/*
if (isnan(X) || isnan(Y) || isnan(Z) || isnan(L) || isnan(U) || isnan(V))  {
    std::cout << X << " " << Y << " " << Z << " " << L << " " << U << " " << V << std::endl;
    std::cout << x << " " << y << std::endl;
    exit(0);
}
*/
			setPixel(x,y,Color(L,U,V));
		}
	}

    ///Computed L component values are in the range [0 to 100].
    ///Computed U component values are in the range [-124 to 220].
    ///Computed V component values are in the range [-140 to 116].

	return;
}

void Image::RGBtoLab()	{
	///Convert to XYZ first
	RGBtoXYZ();

	float REF_X = 95.047f; // Observer= 2°, Illuminant= D65
	float REF_Y = 100.000f;
	float REF_Z = 108.883f;


	for (int y=0;y<height;y++)	{
		for (int x=0;x<width;x++)	{
			float X = getPixel(x,y).r() / REF_X;
			float Y = getPixel(x,y).g() / REF_Y;
			float Z = getPixel(x,y).b() / REF_Z;

			if ( X > 0.008856f ) { X = pow( X , 1.0f/3.0f ); }
			else { X = ( 7.787f * X ) + ( 16.0f/116.0f ); }
			if ( Y > 0.008856f ) { Y = pow( Y , 1.0f/3.0f ); }
			else { Y = ( 7.787f * Y ) + ( 16.0f/116.0f ); }
			if ( Z > 0.008856f ) { Z = pow( Z , 1.0f/3.0f ); }
			else { Z = ( 7.787f * Z ) + ( 16.0f/116.0f ); }

			float l = ( 116.0f * Y ) - 16.0f;
			float a = 500.0f * ( X - Y );
			float b = 200.0f * ( Y - Z );

			setPixel(x,y,Color(l,a,b));
		}
	}
	return;
}

bool onImageBoundary(Image *map, Vector2i const &index, int buffer)  {
	if (outOfBounds(map, index(0)+1, index(1)) ||
		outOfBounds(map, index(0)+1, index(1)+1) ||
		outOfBounds(map, index(0), index(1)+1)) {
		return true;
	}
	return false;
}

Image *Image::subimage(Vector2i const &top_left, Vector2i const &bottom_right)	{
	int sizex = bottom_right(0)-top_left(0);
	int sizey = bottom_right(1)-top_left(1);
	assert(sizex>0 && sizey>0);

	Image *img = new Image(sizex+1,sizey+1,0.0f,0.0f,0.0f,1.0f);
	for (int local_y=0,y=top_left(1);y<=bottom_right(1);y++,local_y++)	{
		for (int local_x=0,x=top_left(0);x<=bottom_right(0);x++,local_x++)	{
			img->setPixel(local_x,local_y, getPixel(x,y));
		}
	}
	return img;
}

float Image::getFocalLengthMillimeters()   {
    return focal_length;
}

bool Image::erase()	{
	///Erase anything previously existing
	for (int y=0;y<height;y++)	{
		delete [] data[y];
	}
	delete [] data;
}

float *Image::getData3Channels()    {
    float *_data = new float[3*width*height];
    int pos = 0;
    for (int y=0;y<height;y++)  {
        for (int x=0;x<width;x++)   {
            Color val = getPixel(x,y);
            _data[pos++] = val.r();
            _data[pos++] = val.g();
            _data[pos++] = val.b();
        }
    }
    return _data;
}

float *Image::getData1Channel(int index)  {
    float *_data = new float[width*height];
    int pos = 0;
    for (int y=0;y<height;y++)  {
        for (int x=0;x<width;x++)   {
            Vector3f val = color2vector3(getPixel(x,y));
            _data[pos++] = val[index];
        }
    }
    return _data;

}

bool Image::flipHorizontal()	{
	for (int y=0;y<height;y++)	{
		for (int x=0;x<width/2;x++)	{
			Color val1 = getPixel(x,y);
			Color val2 = getPixel(width-1-x,y);
			setPixel(x,y,val2);
			setPixel(width-1-x,y,val1);
		}
	}

	return true;
}

bool Image::flipVertical()	{
	for (int y=0;y<height/2;y++)	{
		for (int x=0;x<width;x++)	{
			Color val1 = getPixel(x,y);
			Color val2 = getPixel(x,height-1-y);
			setPixel(x,y,val2);
			setPixel(x,height-1-y,val1);
		}
	}

	return true;
}

#endif
