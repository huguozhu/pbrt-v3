
/*
    pbrt source code is Copyright(c) 1998-2016
                        Matt Pharr, Greg Humphreys, and Wenzel Jakob.

    This file is part of pbrt.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
    IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
    PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */

// core/imageio.cpp*
//
// 本模块提供了图像文件的读写功能，支持以下格式：
// - EXR (OpenEXR): 高动态范围图像格式
// - TGA (Truevision Targa): 8 位真彩色图像格式
// - PNG (Portable Network Graphics): 8 位真彩色图像格式
// - PFM (Portable Float Map): 浮点像素格式
// ReadImage 根据文件扩展名自动选择合适的读取器，
// WriteImage 根据文件扩展名自动选择合适的写入器。
//
#include "imageio.h"
#include "ext/lodepng.h"
#include "ext/targa.h"
#include "fileutil.h"
#include "spectrum.h"

#include <ImfRgba.h>
#include <ImfRgbaFile.h>

namespace pbrt {

// ImageIO Local Declarations
static void WriteImageEXR(const std::string &name, const Float *pixels,
                          int xRes, int yRes, int totalXRes, int totalYRes,
                          int xOffset, int yOffset);
static void WriteImageTGA(const std::string &name, const uint8_t *pixels,
                          int xRes, int yRes, int totalXRes, int totalYRes,
                          int xOffset, int yOffset);
static RGBSpectrum *ReadImageTGA(const std::string &name, int *w, int *h);
static RGBSpectrum *ReadImagePNG(const std::string &name, int *w, int *h);
static bool WriteImagePFM(const std::string &filename, const Float *rgb,
                          int xres, int yres);
static RGBSpectrum *ReadImagePFM(const std::string &filename, int *xres,
                                 int *yres);

// ImageIO Function Definitions

// ReadImage: 根据文件名后缀自动选择合适的格式读取图像
// name - 图像文件路径，支持 .exr/.tga/.png/.pfm 格式
// resolution - 输出参数，返回图像的分辨率（宽度 x 高度）
// 返回值：包含图像像素数据的 RGBSpectrum 数组，失败时返回 nullptr
std::unique_ptr<RGBSpectrum[]> ReadImage(const std::string &name,
                                         Point2i *resolution) {
    if (HasExtension(name, ".exr"))
        return std::unique_ptr<RGBSpectrum[]>(
            ReadImageEXR(name, &resolution->x, &resolution->y));
    else if (HasExtension(name, ".tga"))
        return std::unique_ptr<RGBSpectrum[]>(
            ReadImageTGA(name, &resolution->x, &resolution->y));
    else if (HasExtension(name, ".png"))
        return std::unique_ptr<RGBSpectrum[]>(
            ReadImagePNG(name, &resolution->x, &resolution->y));
    else if (HasExtension(name, ".pfm"))
        return std::unique_ptr<RGBSpectrum[]>(
            ReadImagePFM(name, &resolution->x, &resolution->y));
    Error("Unable to load image stored in format \"%s\" for filename \"%s\".",
          strrchr(name.c_str(), '.') ? (strrchr(name.c_str(), '.') + 1)
                                     : "(unknown)",
          name.c_str());
    return nullptr;
}

// WriteImage: 根据文件名后缀自动选择合适的格式写入图像
// name - 输出图像文件路径，支持 .exr/.pfm/.tga/.png 格式
// rgb - RGB 浮点像素数据（每像素三个 Float 值）
// outputBounds - 输出图像在当前平铺区域中的边界
// totalResolution - 完整图像的总分辨率
void WriteImage(const std::string &name, const Float *rgb,
                const Bounds2i &outputBounds, const Point2i &totalResolution) {
    Vector2i resolution = outputBounds.Diagonal();
    if (HasExtension(name, ".exr")) {
        WriteImageEXR(name, rgb, resolution.x, resolution.y, totalResolution.x,
                      totalResolution.y, outputBounds.pMin.x,
                      outputBounds.pMin.y);
    } else if (HasExtension(name, ".pfm")) {
        WriteImagePFM(name, rgb, resolution.x, resolution.y);
    } else if (HasExtension(name, ".tga") || HasExtension(name, ".png")) {
        // 8-bit formats; apply gamma
        Vector2i resolution = outputBounds.Diagonal();
        std::unique_ptr<uint8_t[]> rgb8(
            new uint8_t[3 * resolution.x * resolution.y]);
        uint8_t *dst = rgb8.get();
        for (int y = 0; y < resolution.y; ++y) {
            for (int x = 0; x < resolution.x; ++x) {
#define TO_BYTE(v) (uint8_t) Clamp(255.f * GammaCorrect(v) + 0.5f, 0.f, 255.f)
                dst[0] = TO_BYTE(rgb[3 * (y * resolution.x + x) + 0]);
                dst[1] = TO_BYTE(rgb[3 * (y * resolution.x + x) + 1]);
                dst[2] = TO_BYTE(rgb[3 * (y * resolution.x + x) + 2]);
#undef TO_BYTE
                dst += 3;
            }
        }

        if (HasExtension(name, ".tga"))
            WriteImageTGA(name, rgb8.get(), resolution.x, resolution.y,
                          totalResolution.x, totalResolution.y,
                          outputBounds.pMin.x, outputBounds.pMin.y);
        else {
            unsigned int error = lodepng_encode24_file(
                name.c_str(), rgb8.get(), resolution.x, resolution.y);
            if (error != 0)
                Error("Error writing PNG \"%s\": %s", name.c_str(),
                      lodepng_error_text(error));
        }
    } else {
        Error("Can't determine image file type from suffix of filename \"%s\"",
              name.c_str());
    }
}

// ReadImageEXR: 读取 OpenEXR 格式的高动态范围图像
// name - EXR 文件路径
// width, height - 输出参数，返回图像尺寸
// dataWindow - 输出参数，返回像素数据窗口（pbrt 使用的非包含性边界约定）
// displayWindow - 输出参数，返回显示窗口
// 返回值：RGBSpectrum 像素数组，失败时返回 NULL
RGBSpectrum *ReadImageEXR(const std::string &name, int *width, int *height,
                          Bounds2i *dataWindow, Bounds2i *displayWindow) {
    using namespace Imf;
    using namespace Imath;
    try {
        RgbaInputFile file(name.c_str());
        Box2i dw = file.dataWindow();

        // OpenEXR uses inclusive pixel bounds; adjust to non-inclusive
        // (the convention pbrt uses) in the values returned.
        // OpenEXR 使用包含性像素边界；调整为 pbrt 使用的非包含性约定
        if (dataWindow)
            *dataWindow = {{dw.min.x, dw.min.y}, {dw.max.x + 1, dw.max.y + 1}};
        if (displayWindow) {
            Box2i dispw = file.displayWindow();
            *displayWindow = {{dispw.min.x, dispw.min.y},
                              {dispw.max.x + 1, dispw.max.y + 1}};
        }
        *width = dw.max.x - dw.min.x + 1;
        *height = dw.max.y - dw.min.y + 1;

        // 分配像素缓冲区并读取图像数据
        std::vector<Rgba> pixels(*width * *height);
        file.setFrameBuffer(&pixels[0] - dw.min.x - dw.min.y * *width, 1,
                            *width);
        file.readPixels(dw.min.y, dw.max.y);

        // 将 OpenEXR 的 RGBA 像素转换为 pbrt 的 RGBSpectrum
        RGBSpectrum *ret = new RGBSpectrum[*width * *height];
        for (int i = 0; i < *width * *height; ++i) {
            Float frgb[3] = {pixels[i].r, pixels[i].g, pixels[i].b};
            ret[i] = RGBSpectrum::FromRGB(frgb);
        }
        LOG(INFO) << StringPrintf("Read EXR image %s (%d x %d)",
                                  name.c_str(), *width, *height);
        return ret;
    } catch (const std::exception &e) {
        Error("Unable to read image file \"%s\": %s", name.c_str(), e.what());
    }

    return NULL;
}

// WriteImageEXR: 以 OpenEXR 格式写入图像数据（内部静态函数）
// name - 输出文件名
// pixels - RGB 浮点像素数据
// xRes, yRes - 当前写入区域的尺寸
// totalXRes, totalYRes - 完整图像的尺寸
// xOffset, yOffset - 当前区域在完整图像中的偏移量
static void WriteImageEXR(const std::string &name, const Float *pixels,
                          int xRes, int yRes, int totalXRes, int totalYRes,
                          int xOffset, int yOffset) {
    using namespace Imf;
    using namespace Imath;

    // 将 Float 像素数据转换为 OpenEXR 的 Rgba 格式
    Rgba *hrgba = new Rgba[xRes * yRes];
    for (int i = 0; i < xRes * yRes; ++i)
        hrgba[i] = Rgba(pixels[3 * i], pixels[3 * i + 1], pixels[3 * i + 2]);

    // OpenEXR uses inclusive pixel bounds.
    Box2i displayWindow(V2i(0, 0), V2i(totalXRes - 1, totalYRes - 1));
    Box2i dataWindow(V2i(xOffset, yOffset),
                     V2i(xOffset + xRes - 1, yOffset + yRes - 1));

    try {
        RgbaOutputFile file(name.c_str(), displayWindow, dataWindow,
                            WRITE_RGB);
        file.setFrameBuffer(hrgba - xOffset - yOffset * xRes, 1, xRes);
        file.writePixels(yRes);
    } catch (const std::exception &exc) {
        Error("Error writing \"%s\": %s", name.c_str(), exc.what());
    }

    delete[] hrgba;
}

// TGA Function Definitions

// WriteImageTGA: 以 TGA 格式写入图像数据
// TGA 格式使用 BGR 像素排列，需要将 RGB 数据重新排列为 BGR
void WriteImageTGA(const std::string &name, const uint8_t *pixels, int xRes,
                   int yRes, int totalXRes, int totalYRes, int xOffset,
                   int yOffset) {
    // Reformat to BGR layout.
    // 将 RGB 重新排列为 TGA 所需的 BGR 格式
    std::unique_ptr<uint8_t[]> outBuf(new uint8_t[3 * xRes * yRes]);
    uint8_t *dst = outBuf.get();
    const uint8_t *src = pixels;
    for (int y = 0; y < yRes; ++y) {
        for (int x = 0; x < xRes; ++x) {
            dst[0] = src[2];
            dst[1] = src[1];
            dst[2] = src[0];
            dst += 3;
            src += 3;
        }
    }

    // 调用 targa 库写入 BGR 格式的 TGA 文件
    tga_result result;
    if ((result = tga_write_bgr(name.c_str(), outBuf.get(), xRes, yRes, 24)) !=
        TGA_NOERR)
        Error("Unable to write output file \"%s\" (%s)",
              name.c_str(), tga_error(result));
}

// ReadImageTGA: 读取 TGA 格式图像文件
// name - TGA 文件路径
// width, height - 输出参数，返回图像尺寸
// 返回值：RGBSpectrum 像素数组，失败时返回 nullptr
static RGBSpectrum *ReadImageTGA(const std::string &name, int *width,
                                 int *height) {
    tga_image img;
    tga_result result;
    if ((result = tga_read(&img, name.c_str())) != TGA_NOERR) {
        Error("Unable to read from TGA file \"%s\" (%s)",
              name.c_str(), tga_error(result));
        return nullptr;
    }

    // 标准化图像方向：统一转为从左到右、从上到下的布局
    if (tga_is_right_to_left(&img)) tga_flip_horiz(&img);
    if (!tga_is_top_to_bottom(&img)) tga_flip_vert(&img);
    if (tga_is_colormapped(&img)) tga_color_unmap(&img);

    *width = img.width;
    *height = img.height;

    // "Unpack" the pixels (origin in the lower left corner).
    // TGA pixels are in BGRA format.
    // 解包像素（原点在左下角），TGA 像素为 BGRA 格式
    RGBSpectrum *ret = new RGBSpectrum[*width * *height];
    RGBSpectrum *dst = ret;
    for (int y = 0; y < *height; y++)
        for (int x = 0; x < *width; x++) {
            uint8_t *src = tga_find_pixel(&img, x, y);
            if (tga_is_mono(&img))
                *dst++ = RGBSpectrum(*src / 255.f);
            else {
                // TGA 是 BGR 格式，转换到 RGB
                Float c[3];
                c[2] = src[0] / 255.f;
                c[1] = src[1] / 255.f;
                c[0] = src[2] / 255.f;
                *dst++ = RGBSpectrum::FromRGB(c);
            }
        }

    tga_free_buffers(&img);
    LOG(INFO) << StringPrintf("Read TGA image %s (%d x %d)",
                              name.c_str(), *width, *height);

    return ret;
}

// ReadImagePNG: 读取 PNG 格式图像文件
// 使用 lodepng 库解码 PNG 文件为 24 位 RGB 数据，然后转换为 RGBSpectrum
static RGBSpectrum *ReadImagePNG(const std::string &name, int *width,
                                 int *height) {
    unsigned char *rgb;
    unsigned w, h;
    unsigned int error = lodepng_decode24_file(&rgb, &w, &h, name.c_str());
    if (error != 0) {
        Error("Error reading PNG \"%s\": %s", name.c_str(),
              lodepng_error_text(error));
        return nullptr;
    }
    *width = w;
    *height = h;

    // 将 8 位 RGB 数据转换为 pbrt 的 RGBSpectrum 格式
    RGBSpectrum *ret = new RGBSpectrum[*width * *height];
    unsigned char *src = rgb;
    for (unsigned int y = 0; y < h; ++y) {
        for (unsigned int x = 0; x < w; ++x, src += 3) {
            Float c[3];
            c[0] = src[0] / 255.f;
            c[1] = src[1] / 255.f;
            c[2] = src[2] / 255.f;
            ret[y * *width + x] = RGBSpectrum::FromRGB(c);
        }
    }

    free(rgb);
    LOG(INFO) << StringPrintf("Read PNG image %s (%d x %d)",
                              name.c_str(), *width, *height);
    return ret;
}

// PFM Function Definitions
/*
 * PFM reader/writer code courtesy Jiawen "Kevin" Chen
 * (http://people.csail.mit.edu/jiawen/)
 */

static PBRT_CONSTEXPR bool hostLittleEndian =
#if defined(__BYTE_ORDER__)
  #if __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
    true
  #elif __BYTE_ORDER__ == __ORDER_BIG_ENDIAN__
    false
  #else
    #error "__BYTE_ORDER__ defined but has unexpected value"
  #endif
#else
  #if defined(__LITTLE_ENDIAN__) || defined(__i386__) || defined(__x86_64__) || \
      defined(_WIN32) || defined(WIN32)
    true
  #elif defined(__BIG_ENDIAN__)
    false
  #elif defined(__sparc) || defined(__sparc__)
    false
  #else
    #error "Can't detect machine endian-ness at compile-time."
  #endif
#endif
    ;

#define BUFFER_SIZE 80

// isWhitespace: 判断字符是否为空白字符（空格、换行、制表符）
static inline int isWhitespace(char c) {
    return c == ' ' || c == '\n' || c == '\t';
}

// readWord: 从文件中读取一个"单词"，即连续的非空白字符序列
// 将读取的字符存入 buffer 并添加 null 终止符。
// 返回读取的字符数（不含终止符和尾随空白），出错时返回 -1
static int readWord(FILE *fp, char *buffer, int bufferLength) {
    int n;
    int c;

    if (bufferLength < 1) return -1;

    n = 0;
    c = fgetc(fp);
    while (c != EOF && !isWhitespace(c) && n < bufferLength) {
        buffer[n] = c;
        ++n;
        c = fgetc(fp);
    }

    if (n < bufferLength) {
        buffer[n] = '\0';
        return n;
    }

    return -1;
}

// ReadImagePFM: 读取 PFM（Portable Float Map）格式图像文件
// 支持单通道（Pf）和三通道（PF）格式，自动处理字节序和缩放因子
static RGBSpectrum *ReadImagePFM(const std::string &filename, int *xres,
                                 int *yres) {
    float *data = nullptr;
    RGBSpectrum *rgb = nullptr;
    char buffer[BUFFER_SIZE];
    unsigned int nFloats;
    int nChannels, width, height;
    float scale;
    bool fileLittleEndian;

    FILE *fp = fopen(filename.c_str(), "rb");
    if (!fp) goto fail;

    // read either "Pf" or "PF"
    if (readWord(fp, buffer, BUFFER_SIZE) == -1) goto fail;

    if (strcmp(buffer, "Pf") == 0)
        nChannels = 1;
    else if (strcmp(buffer, "PF") == 0)
        nChannels = 3;
    else
        goto fail;

    // read the rest of the header
    // read width
    if (readWord(fp, buffer, BUFFER_SIZE) == -1) goto fail;
    width = atoi(buffer);
    *xres = width;

    // read height
    if (readWord(fp, buffer, BUFFER_SIZE) == -1) goto fail;
    height = atoi(buffer);
    *yres = height;

    // read scale
    if (readWord(fp, buffer, BUFFER_SIZE) == -1) goto fail;
    sscanf(buffer, "%f", &scale);

    // read the data
    nFloats = nChannels * width * height;
    data = new float[nFloats];
    // Flip in Y, as P*M has the origin at the lower left.
    for (int y = height - 1; y >= 0; --y) {
        if (fread(&data[y * nChannels * width], sizeof(float),
                  nChannels * width, fp) != nChannels * width)
            goto fail;
    }

    // apply endian conversian and scale if appropriate
    fileLittleEndian = (scale < 0.f);
    if (hostLittleEndian ^ fileLittleEndian) {
        uint8_t bytes[4];
        for (unsigned int i = 0; i < nFloats; ++i) {
            memcpy(bytes, &data[i], 4);
            std::swap(bytes[0], bytes[3]);
            std::swap(bytes[1], bytes[2]);
            memcpy(&data[i], bytes, 4);
        }
    }
    if (std::abs(scale) != 1.f)
        for (unsigned int i = 0; i < nFloats; ++i) data[i] *= std::abs(scale);

    // create RGBs...
    rgb = new RGBSpectrum[width * height];
    if (nChannels == 1) {
        for (int i = 0; i < width * height; ++i) rgb[i] = RGBSpectrum(data[i]);
    } else {
        for (int i = 0; i < width * height; ++i) {
            Float frgb[3] = {data[3 * i], data[3 * i + 1], data[3 * i + 2]};
            rgb[i] = RGBSpectrum::FromRGB(frgb);
        }
    }

    delete[] data;
    fclose(fp);
    LOG(INFO) << StringPrintf("Read PFM image %s (%d x %d)",
                              filename.c_str(), *xres, *yres);
    return rgb;

fail:
    Error("Error reading PFM file \"%s\"", filename.c_str());
    if (fp) fclose(fp);
    delete[] data;
    delete[] rgb;
    return nullptr;
}

// WriteImagePFM: 以 PFM 格式写入三通道（PF）浮点图像
// PFM 格式将像素从底部左到顶部右逐行存储，使用 32 位浮点数
static bool WriteImagePFM(const std::string &filename, const Float *rgb,
                          int width, int height) {
    FILE *fp;
    float scale;

    fp = fopen(filename.c_str(), "wb");
    if (!fp) {
        Error("Unable to open output PFM file \"%s\"", filename.c_str());
        return false;
    }

    std::unique_ptr<float[]> scanline(new float[3 * width]);

    // only write 3 channel PFMs here...
    // 仅写入三通道 PFM 格式（标记为 "PF"）
    if (fprintf(fp, "PF\n") < 0) goto fail;

    // write the width and height, which must be positive
    // 写入宽度和高度
    if (fprintf(fp, "%d %d\n", width, height) < 0) goto fail;

    // write the scale, which encodes endianness
    // 写入缩放因子，其正负号编码了字节序（负值表示小端）
    scale = hostLittleEndian ? -1.f : 1.f;
    if (fprintf(fp, "%f\n", scale) < 0) goto fail;

    // write the data from bottom left to upper right as specified by
    // http://netpbm.sourceforge.net/doc/pfm.html
    // The raster is a sequence of pixels, packed one after another, with no
    // delimiters of any kind. They are grouped by row, with the pixels in each
    // row ordered left to right and the rows ordered bottom to top.
    // 按 PFM 规范从下到上写入行数据
    for (int y = height - 1; y >= 0; y--) {
        // in case Float is 'double', copy into a staging buffer that's
        // definitely a 32-bit float...
        // 如果 Float 是 double 类型，先复制到 32 位 float 缓冲区
        for (int x = 0; x < 3 * width; ++x)
            scanline[x] = rgb[y * width * 3 + x];
        if (fwrite(&scanline[0], sizeof(float), width * 3, fp) <
            (size_t)(width * 3))
            goto fail;
    }

    fclose(fp);
    return true;

fail:
    Error("Error writing PFM file \"%s\"", filename.c_str());
    fclose(fp);
    return false;
}

}  // namespace pbrt
