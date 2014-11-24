/**
* Copyright (C) 2013 Raúl Moreno Galdón
*
* This program is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "data_collector.h"

#ifdef CHECK_DUPLICATES
	static char *ult_seq = NULL;
	static U_CYCLES l_ult_seq = 0;
	static uint64_t pos_ult_seq;
#endif

/**
 * Get recalibration data from alignment.
 */
ERROR_CODE
recal_get_data_from_bam_alignment(const bam1_t* read, const genome_t* ref, recal_info_t* output_data, recal_data_collect_env_t *collect_env)
{
	assert(read);
	assert(ref);
	assert(output_data);
	assert(collect_env);

	char *ref_seq;
	bam1_t *prev_read;

	int64_t init_pos, end_pos;
	size_t init_pos_ref, end_pos_ref;
	char *comp_res;
	char *comp_mask;
	char *dinucs;
	uint32_t flag;
	unsigned int chr;

	//Enviroment
	char *bam_seq;
	char *bam_quals;
	U_CYCLES bam_seq_l;

	//uint32_t *read_cigar;
	size_t read_disp_clip;
	size_t ref_disp;

	char *read_seq_ref;
	size_t read_seq_ref_l;
	uint32_t misses;
	size_t read_l;
	uint32_t tmp_cigar[30];

	unsigned int i;

	//This read has been filtered?
	if(read->core.flag & 2048)
		return NO_ERROR;

	//SET VARS
	{
		bam_seq = collect_env->bam_seq;
		bam_quals = collect_env->bam_quals;
		bam_seq_l = 0;
		prev_read = collect_env->prev_read;
	}

	//Check duplicate
	if(prev_read)
	{
		if(prev_read->core.tid == read->core.tid && prev_read->core.pos == read->core.pos)
			return NO_ERROR;
	}
	collect_env->prev_read = read;

	//Check duplicate
	if(prev_read)
	{
		if(prev_read->core.tid == read->core.tid && prev_read->core.pos == read->core.pos)
			return NO_ERROR;
	}
	collect_env->prev_read = read;

	//Get sequence
	new_sequence_from_bam_ref((bam1_t *)read, bam_seq, read->core.l_qseq + 1);

	//Get quals
	new_quality_from_bam_ref((bam1_t *)read, 0, bam_quals, read->core.l_qseq + 1);

	//LOG_WARN("========\n");
	//LOG_WARN_F("SEQ:%3d - %s\n", read->core.l_qseq, bam_seq);

	//Hardclip
	cigar32_hardclip_softclips(bam1_cigar(read), read->core.n_cigar, bam_seq, bam_quals, read->core.l_qseq, tmp_cigar, bam_seq, bam_quals, NULL);

	//Get cycles and positions
	//cycles = alig->core.l_qseq;
	//bam_seq_l = aux_res_seq_l;
	//bam_seq_l = read->core.l_qseq;
	bam_seq_l = strlen(bam_seq);

	//LOG_WARN_F("HCP:%3d - %s\n", bam_seq_l, bam_seq);

	if(bam_seq_l == 0)
	{
		LOG_WARN_F("Alignment with sequence length zero: %s\n", bam1_qname(read));
		return NO_ERROR;
	}
	ref_disp = 100;
	init_pos = read->core.pos - ref_disp;
	if(init_pos < 0)
	{
		ref_disp = read->core.pos - abs(init_pos);
		init_pos = 1;
	}

	end_pos = init_pos + (bam_seq_l  * 4) + 100;

	init_pos_ref = init_pos + RECAL_REFERENCE_CORRECTION_OFFSET;
	end_pos_ref = end_pos + RECAL_REFERENCE_CORRECTION_OFFSET;

	//Duplicates check
	#ifdef CHECK_DUPLICATES
	{
		if(ult_seq == NULL)
			ult_seq = (char *)malloc(sizeof(char) * output_data->num_cycles);

		if(l_ult_seq)
		{
			if(pos_ult_seq == read->core.pos && l_ult_seq == bam_seq_l && strcmp(bam_seq, ult_seq) == 0)
			{
				//LOG_WARN_F("\nDUPLICATE POS: %d CYCLES: %d\n\tSEQ:  %s\n\tLAST: %s", init_pos, bam_seq_l, bam_seq, ult_seq);
				duplicated++;
				return NO_ERROR;
			}
		}
	}
	#endif

	//Obtain reference for this 100 nucleotides
	flag = (uint32_t) read->core.flag;
	chr = read->core.tid;
	if((flag & BAM_FUNMAP) || init_pos_ref == 0 || end_pos_ref == 0 || ref->chr_size[chr] == 0)
	{
		//Read is unmapped
		//LOG_WARN_F("Alignment is unmapped %d:%d-%d, %s\n", read->core.tid + 1, init_pos_ref + 1, end_pos_ref + 1, bam1_qname(read));
		return NO_ERROR;
	}

	ref_seq = (char *)malloc(((end_pos_ref - init_pos_ref) + 2) * sizeof(char));
	genome_read_sequence_by_chr_index(ref_seq, (flag & BAM_FREVERSE) ? 1 : 0, chr, &init_pos_ref, &end_pos_ref, (genome_t *)ref);

	if((end_pos_ref - init_pos_ref) == 0)
	{
		LOG_WARN_F("Cannot obtain reference for region %d:%lu-%lu, %s\n", read->core.tid + 1, init_pos_ref + 1, end_pos_ref + 1, bam1_qname(read));
		free(ref_seq);
		return NO_ERROR;
	}

	//Allocations
	comp_res = (char *)malloc((read->core.l_qseq + 1) * sizeof(char));
	comp_mask = (char *)malloc((read->core.l_qseq + 1) * sizeof(char));
	dinucs = (char *)malloc((read->core.l_qseq + 1) * sizeof(char));
	memset(comp_res, 0, (read->core.l_qseq + 1) * sizeof(char));
	memset(comp_mask, 0, (read->core.l_qseq + 1) * sizeof(char));

	//Get initial clip displacement
	//cigar32_count_clip_displacement(tmp_cigar, read->core.n_cigar, &read_disp_clip);
	cigar32_count_nucleotides_not_clip(tmp_cigar, read->core.n_cigar, &read_l);

	//Create sequence to compare with reference
	read_seq_ref = (char *) malloc((read->core.l_qseq + 1) * sizeof(char));
	cigar32_create_ref(tmp_cigar, read->core.n_cigar,
			ref_seq + ref_disp, strlen(ref_seq) - ref_disp,
			bam_seq, strlen(bam_seq),
			read_seq_ref, &read_seq_ref_l, comp_mask);
	//LOG_WARN_F("REF:%3d - %s\n", strlen(read_seq_ref), read_seq_ref);

	//Correct comparation?
	if(read_seq_ref_l != strlen(bam_seq))
	{
		LOG_WARN_F("Read-Ref: %d, Read: %d, Not enough reference sequence length\n", read_seq_ref_l, strlen(bam_seq));
	}

	//Get raw score with reference
	nucleotide_compare(read_seq_ref, bam_seq, read_seq_ref_l, comp_res, &misses);

	//Avoid N nucleotides
	#ifdef NOT_COUNT_NUCLEOTIDE_N
	for(i = 0; i < read_seq_ref_l; i++)
	{
		if(read_seq_ref[i] == 'N')
		{
			misses--;
			comp_mask[i] = 0;
		}
	}
	#endif

	/*LOG_WARN_F("CMP:%3d - ", read_seq_ref_l);
	for(i=0;i < read_seq_ref_l; i++)
		fprintf(log_file, "%c", comp_res[i] == 0 ? '0' : '1');
	fprintf(log_file, "\n");*/

	//Dinucs
	for(i = 0; i < bam_seq_l; i++)
	{
		if(i > 0)
		{
		  recal_get_dinuc(bam_seq[i-1], bam_seq[i], (U_DINUC *)&dinucs[i]);
		}
		else
		{
			dinucs[i] = d_X;
		}
	}

	//Add data
	recal_add_base_v(output_data, bam_seq, bam_quals, 0, bam_seq_l, dinucs, comp_res, comp_mask);
	/*LOG_WARN_F("MSK:%3d - ", bam_seq_l);
		for(i=0;i < bam_seq_l; i++)
			fprintf(log_file, "%c", comp_mask[i] == 0 ? '0' : '1');
		fprintf(log_file, "\n");
	LOG_WARN_F("MSS:%3d - ", bam_seq_l);
			for(i=0;i < bam_seq_l; i++)
				fprintf(log_file, "%c", comp_res[i] == 0 ? '0' : '1');
			fprintf(log_file, "\n");*/

	//Set last sequence for duplicates
	#ifdef CHECK_DUPLICATES
	{
		strcpy(ult_seq, bam_seq);
		l_ult_seq = read->core.l_qseq;
		pos_ult_seq = read->core.pos;
	}
	#endif

	//Free resources
	{
		free(ref_seq);
		free(comp_res);
		free(comp_mask);
		free(dinucs);
		free(read_seq_ref);
	}

	return NO_ERROR;
}
