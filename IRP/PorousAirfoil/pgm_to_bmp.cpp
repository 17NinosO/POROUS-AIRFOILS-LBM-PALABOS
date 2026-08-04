// pgm_to_bmp.cpp
#include <fstream>
#include <vector>
#include <cstdint>
#include <iostream>
#include <stdexcept>

void writeLE32(std::ofstream& out, uint32_t v) { out.write(reinterpret_cast<char*>(&v), 4); }
void writeLE16(std::ofstream& out, uint16_t v) { out.write(reinterpret_cast<char*>(&v), 2); }

int main(int argc, char** argv) {
    if (argc != 3) {
        std::cerr << "usage: pgm_to_bmp input.pgm output.bmp\n";
        return 1;
    }

    std::ifstream in(argv[1], std::ios::binary);
    if (!in) throw std::runtime_error("cannot open input PGM");

    std::string magic;
    int width, height, maxval;
    in >> magic >> width >> height >> maxval;
    in.get();  // consume the single whitespace char before binary data

    if (magic != "P5") throw std::runtime_error("only binary P5 PGM supported");

    std::vector<unsigned char> pixels(size_t(width) * height);
    in.read(reinterpret_cast<char*>(pixels.data()), pixels.size());

    int row_padded = (width + 3) & ~3;  // BMP rows are padded to a multiple of 4 bytes
    uint32_t pixel_data_offset = 14 + 40 + 256 * 4;
    uint32_t file_size = pixel_data_offset + row_padded * height;

    std::ofstream out(argv[2], std::ios::binary);

    // --- file header ---
    out.write("BM", 2);
    writeLE32(out, file_size);
    writeLE32(out, 0);
    writeLE32(out, pixel_data_offset);

    // --- DIB header (BITMAPINFOHEADER) ---
    writeLE32(out, 40);
    writeLE32(out, uint32_t(width));
    writeLE32(out, uint32_t(height));
    writeLE16(out, 1);
    writeLE16(out, 8);
    writeLE32(out, 0);
    writeLE32(out, row_padded * height);
    writeLE32(out, 0);
    writeLE32(out, 0);
    writeLE32(out, 256);
    writeLE32(out, 0);

    // --- greyscale palette (256 shades) ---
    for (int i = 0; i < 256; ++i) {
        out.put(char(i)); out.put(char(i)); out.put(char(i)); out.put(0);
    }

    // --- pixel data, bottom row first (BMP convention) ---
    std::vector<char> pad(row_padded - width, 0);
    for (int row = height - 1; row >= 0; --row) {
        out.write(reinterpret_cast<char*>(&pixels[size_t(row) * width]), width);
        if (!pad.empty()) out.write(pad.data(), pad.size());
    }

    std::cout << "wrote " << argv[2] << " (" << width << "x" << height << ")\n";
    return 0;
}