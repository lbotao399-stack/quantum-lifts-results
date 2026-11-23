// Exact fixed-content necklaces and canonical 14-letter codes. No linear
// combination, finite-colour relation or change of coordinates is introduced.
#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <exception>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
thread_local std::string error_message;

struct Catalogue {
    int length;
    std::vector<int> symbols, remaining, word;
    std::vector<std::uint64_t> result;

    Catalogue(int alphabet, const int* counts) : length(0) {
        if (alphabet < 1 || alphabet > 4)
            throw std::invalid_argument("alphabet must contain 1 through 4 symbols");
        for (int i = 0; i < alphabet; ++i) {
            if (counts[i] < 0 || counts[i] > 31 || length > 31 - counts[i])
                throw std::invalid_argument("content outside supported trace length");
            length += counts[i];
            if (counts[i]) { symbols.push_back(i); remaining.push_back(counts[i]); }
        }
        if (length < 2 || length > 31)
            throw std::invalid_argument("trace length must be 2 through 31");
        word.assign(length + 1, 0);
        --remaining[0];
    }

    void emit() {
        std::uint64_t forward = 0, reverse = 0;
        for (int i = 1; i <= length; ++i) {
            forward = (forward << 2) | symbols[word[i]];
            reverse = (reverse << 2) | symbols[word[length + 1 - i]];
        }
        const auto mask = (std::uint64_t(1) << (2 * length)) - 1;
        auto minimum_reverse = reverse;
        for (int i = 1; i < length; ++i) {
            reverse = ((reverse << 2) & mask) | (reverse >> (2 * length - 2));
            minimum_reverse = std::min(minimum_reverse, reverse);
        }
        if (forward > minimum_reverse || (length % 2 && forward == minimum_reverse))
            return;
        result.push_back(forward);
    }

    void visit(int position, int period) {
        if (position > length) {
            if (length % period == 0) emit();
            return;
        }
        const int lower = word[position - period];
        for (int letter = lower; letter < int(symbols.size()); ++letter) {
            if (remaining[letter]) {
                word[position] = letter;
                --remaining[letter];
                visit(position + 1, letter == lower ? period : position);
                ++remaining[letter];
            }
        }
    }
};
}

extern "C" {
const char* dihedral_error() { return error_message.c_str(); }

std::uint64_t* dihedral_catalogue(int alphabet, const int* counts, std::size_t* size) {
    try {
        error_message.clear();
        Catalogue catalogue(alphabet, counts);
        catalogue.visit(2, 1);
        std::sort(catalogue.result.begin(), catalogue.result.end());
        *size = catalogue.result.size();
        auto* output = static_cast<std::uint64_t*>(std::malloc(sizeof(std::uint64_t) * std::max<std::size_t>(1, *size)));
        if (!output) throw std::bad_alloc();
        std::copy(catalogue.result.begin(), catalogue.result.end(), output);
        return output;
    } catch (const std::exception& error) {
        error_message = error.what();
        *size = 0;
        return nullptr;
    }
}

void dihedral_free(void* memory) { std::free(memory); }

void dihedral_canonical14(const std::uint32_t* input, std::uint32_t* output,
                          std::size_t size) {
    constexpr std::uint32_t mask = (1u << 28) - 1;
    for (std::size_t j = 0; j < size; ++j) {
        auto forward = input[j];
        std::uint32_t reverse = 0;
        for (int i = 0; i < 14; ++i)
            reverse = (reverse << 2) | ((forward >> (2 * i)) & 3);
        auto smallest = std::min(forward, reverse);
        for (int i = 1; i < 14; ++i) {
            forward = ((forward << 2) & mask) | (forward >> 26);
            reverse = ((reverse << 2) & mask) | (reverse >> 26);
            smallest = std::min(smallest, std::min(forward, reverse));
        }
        output[j] = smallest;
    }
}
}
