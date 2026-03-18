//
// Copyright (C) 2025 John Auld
//

// C++23
//
// Compile with:
//
// g++ --std=c++23 -Wall -Wextra -O3 -march=native -fopenmp probeFuzz.cpp
//

#include <vector>
#include <iostream>
#include <mutex>
#include <sstream>

namespace {

constexpr size_t	KMER_LEN	{50};
constexpr size_t	LO_BASES	{32};					// fits in one uint64_t
constexpr size_t	HI_BASES	{KMER_LEN - LO_BASES};	// 18
constexpr int		MAXD		{10};

// Encode A,C,G,T into 2 bits.
inline uint64_t enc(char c) {
	switch (c) {
	case 'A': case 'a' : return 0;
	case 'C': case 'c' : return 1;
	case 'G': case 'g' : return 2;
	case 'T': case 't' : return 3;
	default: throw std::domain_error(std::string{"Invalid base: "} + c);
	}
}

// Stores a 50-mer as two uint64_t words (32 bases in lo_, the rest in hi_)
struct Packed {
	uint64_t lo_;
	uint64_t hi_;
};

[[nodiscard]]
inline Packed pack50(std::string_view s) {
	if (s.size() < KMER_LEN)
		throw std::domain_error("Expected 50-mer, got: " + std::string{s});

	Packed p {0, 0};

	for ( size_t i = 0; i < LO_BASES; ++i )
		p.lo_ |= enc(s[i]) << (i * 2);

	for ( size_t i = 0; i < HI_BASES; ++i )
		p.hi_ |= enc(s[i + LO_BASES]) << (i * 2);

	return p;
}

[[nodiscard]]
inline std::string unpack50(Packed packed) {	// by value, since we munge it
	static constexpr char bases[] { 'a', 'c', 'g', 't' };

	std::string ret;
	ret.reserve(KMER_LEN);

	size_t count {LO_BASES};
	while ( count-- ) {
		ret.push_back(bases[packed.lo_ & 0x3]);
		packed.lo_ >>= 2;
	}	

	count = HI_BASES;
	while ( count-- ) {
		ret.push_back(bases[packed.hi_ & 0x3]);
		packed.hi_ >>= 2;
	}

	return ret;
}

// Compute Hamming distance of two packed 50-mers
[[nodiscard]]
inline int hamming(const Packed &a, const Packed &b) {
	// Get the diffs between a and b (as 2-bit groups since each base is 2 bits)
	uint64_t lo {a.lo_ ^ b.lo_};
	uint64_t hi {a.hi_ ^ b.hi_};

	// Collapse 2-bit groups into 1-bit mismatch flags
	lo = (lo | (lo >> 1)) & 0x5555'5555'5555'5555ULL;
	hi = (hi | (hi >> 1)) & 0x5555'5555'5555'5555ULL;

	// Count the 1 bits in each group and return the sum
	return __builtin_popcountll(lo) + __builtin_popcountll(hi);
}

} // anonymous namespace

int main() {
	std::ios::sync_with_stdio(false);
	std::cin.tie(nullptr);

	std::vector<std::string>	id;
	std::vector<Packed>			seq;
	std::string					line;

	while ( std::cin >> line ) {					// Read until EOF
		auto comma_pos {line.find(',')};

		if ( comma_pos == std::string::npos ) continue;	// Skip malformed lines

		id.push_back(line.substr(0, comma_pos));
		seq.push_back(pack50(line.substr(comma_pos + 1)));
	}

	const auto	num_probes	{seq.size()};

	std::mutex	cout_mtx;

	#pragma omp parallel
	{
		std::ostringstream out;

		#pragma omp for schedule(dynamic, 64)
		for ( size_t i = 0; i < num_probes; ++i )
			for ( size_t j = i+1; j < num_probes; ++j )
				if ( auto d = hamming(seq[i], seq[j]); d <= MAXD )
					out << d << ',' << id[i] << ',' << id[j] << ','
						<< unpack50(seq[i]) << ',' << unpack50(seq[j]) << '\n';

		std::scoped_lock lock {cout_mtx};
		std::cout << out.str();
	}	// omp parallel
}
