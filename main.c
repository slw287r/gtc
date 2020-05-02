#include <stdio.h>
#include <stdint.h>
#include <inttypes.h>
#include <assert.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/ioctl.h>

#include "htslib/ketopt.h"
#include "htslib/kseq.h"
#include "htslib/khashl.h"
#include "htslib/cgranges.h"
#include "htslib/faidx.h"
#include "2bit.h"
#include "zlib.h"

KSEQ_INIT(gzFile, gzread)
KHASHL_MAP_INIT(, kc_c1_t, kc_c1, uint64_t, uint32_t, kh_hash_uint64, kh_eq_generic)

#define PRIMARY 24 // primary chromosomes
#define max(a, b) (((a) < (b)) ? (b) : (a))
#define min(a, b) (((a) < (b)) ? (a) : (b))
#define basename(str) (strrchr(str, '/') ? strrchr(str, '/') + 1 : str)

const uint8_t kmertochar[5] = { 'A', 'C', 'G', 'T', 'N' };

/**
 * @brief translate ACGT to 0123
 */
const unsigned char seq_nt4_table[256] = {
	0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

/**
 * @brief trinucleotides in output
 */
char sig3[][4] = {
	"ACA", "ACC", "ACG", "ACT",
	"CCA", "CCC", "CCG", "CCT",
	"GCA", "GCC", "GCG", "GCT",
	"TCA", "TCC", "TCG", "TCT",
	"ATA", "ATC", "ATG", "ATT",
	"CTA", "CTC", "CTG", "CTT",
	"GTA", "GTC", "GTG", "GTT",
	"TTA", "TTC", "TTG", "TTT"
};

/**
 * @brief trinucleotides used by kmer_cnt
 */
char key3[][4] = {
	"ACA", "ACC", "ACG", "ACT",
	"CCA", "CCC", "CCG", "AGG",
	"GCA", "GCC", "CGC", "AGC",
	"TCA", "GGA", "CGA", "AGA",
	"ATA", "ATC", "ATG", "AAT",
	"CTA", "CTC", "CAG", "AAG",
	"GTA", "GAC", "CAC", "AAC",
	"TAA", "GAA", "CAA", "AAA"
};

static void usage(char *str);

/**
 * @brief seq to uint64_t
 */
static uint64_t s64u(const char* s)
{
	int i, l;
	uint64_t x[2], mask = (1ULL<<6) - 1, shift = 4;
	for (i = 0, x[0] = x[1] = 0; i < 3; ++i)
	{
		int c = seq_nt4_table[(uint8_t)s[i]];
		x[0] = (x[0] << 2 | c) & mask;                  // forward strand
		x[1] = x[1] >> 2 | (uint64_t)(3 - c) << shift;  // reverse strand
	}
	return x[0] < x[1]? x[0] : x[1];
}

// convert uint64_t to seq
static void u64s(uint64_t x, char* s)
{
	size_t i = 0;
	for (i = 0; i < 3; ++i)
	{
		s[2 - i] = "ACGT"[(uint8_t)x & 0x3];
		x >>= 2;
	}
	s[i] = '\0';
}

// convert uint64_t to reverse complementary seq
static void u64rc(uint64_t x, char* s)
{
	size_t i = 0;
	for (i = 2; i > 0; --i)
	{
		s[2 - i] = "ACGT"[3 - (uint8_t)x & 0x3];
		x >>= 2;
	}
	s[2] = "ACGT"[3 - (uint8_t)x & 0x3];
	s[3] = '\0';
}

// insert k-mers to hash table
static void count_seq(kc_c1_t *h, int k, int len, char *seq)
{
	int i, l;
	uint64_t x[2], mask = (1ULL<<k*2) - 1, shift = (k - 1) * 2;
	for (i = l = 0, x[0] = x[1] = 0; i < len; ++i)
	{
		int absent, c = seq_nt4_table[(uint8_t)seq[i]];
		if (c < 4) // not an "N" base
		{
			x[0] = (x[0] << 2 | c) & mask;                  // forward strand
			x[1] = x[1] >> 2 | (uint64_t)(3 - c) << shift;  // reverse strand
			if (++l >= k) // we find a k-mer
			{
				khint_t itr;
				uint64_t y = x[0] < x[1]? x[0] : x[1];
				itr = kc_c1_put(h, y, &absent); // only add one strand!
				if (absent) kh_val(h, itr) = 0;
				++kh_val(h, itr);
			}
		} else l = 0, x[0] = x[1] = 0; // if there is an "N", restart
	}
}

// get kmer count from hash table
static void get_count(const kc_c1_t *h, const char *s, int *c)
{
	khint_t k;
	for (k = 0; k < kh_end(h); ++k)
	{
		if (kh_exist(h, k))
		{
			if (s64u(s) == kh_key(h, k))
			{
				*c = kh_val(h, k);
				break;
			}
		}
	}
}

/*
 * https://github.com/parklab/SigMA/blob/master/R/get_trinuc_norm.R
 *
 * The get_trinuc_norm() can be used to determine the
 * normalization to be used to match the frequency in the
 * sequencing platform to the frequency in the whole genomes
 *
 * @param bed_file the path to the bed file that defines the
 * library of the sequencing platform
 * @param do_MC_sampling if set to TRUE, to speed up
 * the calculation the regions are sampled randomly rather
 * than providing a full count
 */
static void get_trinuc_norm(const kc_c1_t *h, const kc_c1_t *h_wg, FILE *fp)
{
	int i, j;
	/*
	 * counts_expanded <- c(rep(counts[1:16], 3), rep(counts[17:32], 3))
	 * norm <- counts_expanded/counts_trinuc_genome
	 * return(norm=100*norm/sum(norm))
	 */
	double norm[32] = {0}, sum = 0;
	int c[32] = {0}, c_wg[32] = {0};
	for (i = 0; i < 32; ++i)
	{
		get_count(h, key3[i], c + i);
		get_count(h_wg, key3[i], c_wg + i);
		norm[i] = (double)c[i] / (c_wg[i] ? c_wg[i] : 1);
		sum += norm[i] * 3;
	}
	sum = sum ? sum : 1;
	for (j = 0; j < 3; ++j)
		for (i = 0; i < 16; ++i)
			fprintf(fp, "%s\t%f\n", sig3[i], norm[i] / sum * 100);
	for (j = 0; j < 3; ++j)
		for (i = 16; i < 32; ++i)
			fprintf(fp, "%s\t%f\n", sig3[i], norm[i] / sum * 100);
}

// output in the order of get_trinuc_norm.R
static void print_sig3(const kc_c1_t *h, FILE *fp)
{
	int i;
	for (i = 0; i < 32; ++i)
	{
		int c = 0;
		get_count(h, key3[i], &c);
		fprintf(fp, "%s\t%d\n", sig3[i], c);
	}
}

static char *parse_bed3(char *s, int32_t *st_, int32_t *en_)
{
	char *p, *q, *ctg = 0;
	int32_t i, st = -1, en = -1;
	for (i = 0, p = q = s;; ++q)
	{
		if (*q == '\t' || *q == '\0')
		{
			int c = *q;
			*q = 0;
			if (i == 0) ctg = p;
			else if (i == 1) st = atol(p);
			else if (i == 2) en = atol(p);
			++i, p = q + 1;
			if (i == 3 || c == '\0') break;
		}
	}
	*st_ = st, *en_ = en;
	return i >= 3? ctg : 0;
}

// load panel region bed file to hash
static cgranges_t *read_bed(const char *fn)
{
	gzFile fp;
	kstream_t *ks;
	cgranges_t *cr;
	kstring_t str = {0,0,0};
	int64_t k = 0;
	if ((fp = gzopen(fn, "r")) == 0) return 0;
	ks = ks_init(fp);
	cr = cr_init();
	while (ks_getuntil(ks, KS_SEP_LINE, &str, 0) >= 0)
	{
		char *ctg;
		int32_t st, en;
		ctg = parse_bed3(str.s, &st, &en);
		if (ctg) cr_add(cr, ctg, st, en, k++);
	}
	if (k > INT32_MAX)
		fprintf(stderr, "WARNING: more than %d records; some functionality may not work!\n", INT32_MAX);
	free(str.s);
	ks_destroy(ks);
	gzclose(fp);
	return cr;
}

void context(char *ref, const char *bed, const char *out)
{
	size_t i;
	char *chr;
	int beg, end, ret;
	int is2bit = is_twobit(ref);
	cgranges_t *cr = read_bed(bed); assert(cr);
	if (!cr_is_sorted(cr)) cr_sort(cr);
	//cr_merge_pre_index(cr);
	FILE *fp = out ? fopen(out, "w") : stdout;
	//fputs("#kmer\trevcom\tcount\n", fp);
	faidx_t *fai = is2bit ? NULL : fai_load(ref);
	TwoBit *fb = is2bit ? twobitOpen(ref, 0) : NULL;
	kc_c1_t *h = kc_c1_init();
	// iterate over the merged regions
	for (i = 0; i < cr->n_r; ++i)
	{
		chr = cr->ctg[cr->r[i].x >> 32].name;
		const int len = is2bit ? twobitChromLen(fb, chr) : faidx_seq_len(fai, chr);
		beg = max(0, (int32_t)cr->r[i].x - 2);
		end = min(len, cr->r[i].y + 1);
		char *seq = is2bit ? twobitSequence(fb, chr, beg, end) : faidx_fetch_seq(fai, chr, beg, end - 1, &ret);
		count_seq(h, 3, end - beg, seq);
		free(seq);
	}
	kc_c1_t *h_wg = kc_c1_init();
	//int n = is2bit ? fb->hdr->nChroms : faidx_nseq(fai);
	for (i = 0; i < PRIMARY; ++i) // primary chromosomes only
	{
		const char *chr = is2bit ? fb->cl->chrom[i] : faidx_iseq(fai, i);
		const int len = is2bit ? twobitChromLen(fb, chr) : faidx_seq_len(fai, chr);
		char *seq = is2bit ? twobitSequence(fb, chr, 0, len) : faidx_fetch_seq(fai, chr, 0, len - 1, &ret);
		count_seq(h_wg, 3, len, seq);
		free(seq);
	}
	if (is2bit) twobitClose(fb);
	else fai_destroy(fai);
	cr_destroy(cr);
	// output
	get_trinuc_norm(h, h_wg, fp);
	kc_c1_destroy(h);
	kc_c1_destroy(h_wg);
	fclose(fp);
}

int main(int argc, char *argv[])
{
	int c;
	char *b = 0, *o = 0, *r = 0;
	ketopt_t opt = KETOPT_INIT;
	while ((c = ketopt(&opt, argc, argv, 1, "b:o:", 0)) >= 0)
	{
		if (c == 'b') b = opt.arg;
		if (c == 'o') o = opt.arg;
	}
	if (argc - opt.ind < 1)
	{
		puts("[ERROR] reference file (fa or 2bit) required!");
		usage(argv[0]);
	}
	if (!b || (access(b, R_OK) == -1))
	{
		fputs("[ERROR] required bed file is not specified or cannot be accessed!\n", stderr);
		exit(1);
	}
	r = argv[opt.ind];
	if ((access(r, R_OK) == -1))
	{
		fputs("[ERROR] reference file specified cannot be accessed!\n", stderr);
		exit(1);
	}
	context(r, b, o);
	return 0;
}

static void usage(char *str)
{
	putchar('\n');
	puts("Calculate weight (norm96) for gene panel bed for SigMA");
	putchar('\n');
	fprintf(stdout, "Usage: \e[1;31m%s\e[0;0m [options] <ref.fa | ref.2bit>\n", basename(str));
	putchar('\n');
	puts("  -b gene panel bed file");
	puts("  -o norm96 relative to genome [stdout]");
	puts("\nNotes: For fasta (fa) reference, fai index is required.");
	puts("       You can also make use of BSgenome's 2bit file in");
	puts("       R> \e[3msystem.file(\"extdata/single_sequences.2bit\",");
	puts("                 package=\"BSgenome.Hsapiens.UCSC.hg19\")\e[0m");
	putchar('\n');
	exit(EXIT_FAILURE);
}
