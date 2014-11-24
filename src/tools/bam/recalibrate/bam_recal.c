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

#include "bam_recal.h"

/**
 * Recalibrate alignment and store in file.
 */
ERROR_CODE
recal_recalibrate_alignment(bam1_t* alig, const recal_info_t *bam_info, recal_recalibration_env_t *recalibration_env)
{
	U_QUALS qual_index;
	unsigned int matrix_index;
	unsigned int i;
	U_DINUC dinuc;

	//Sequence
	char *bam_seq;
	char *bam_quals;
	char *res_quals;
	U_CYCLES bam_seq_l;
	U_CYCLES bam_seq_max_l;

	//Recalibration
	double delta_r, delta_rc, delta_rd;

	//CHECK ARGUMENTS (Assuming this function is called always from recal_recalibrate_batch)
	{
		//Check nulls
		ASSERT(alig);
		ASSERT(bam_info);
		ASSERT(recalibration_env);
	}

	//SET VARS
	{
		bam_quals = recalibration_env->bam_quals;
		bam_seq_l = 0;
		bam_seq_max_l = recalibration_env->bam_seq_max_l;
	}

	//Get sequence length
	bam_seq_l = alig->core.l_qseq;

	//Sequence length check
	ASSERT(alig->core.l_qseq <= bam_seq_max_l);
	ASSERT(bam_seq_l != 0);

	//Get sequence
	bam_seq = new_sequence_from_bam((bam1_t *)alig);

	//Get quals
	new_quality_from_bam_ref((bam1_t *)alig, 0, bam_quals, bam_seq_max_l);

	//Allocate for result
	res_quals = (char *)malloc(bam_seq_l * sizeof(char));
	memcpy(res_quals, bam_quals, bam_seq_l * sizeof(char));

	//Iterates nucleotides in this read
	dinuc = 0;
	for(i = 0; i < bam_seq_l; i++)
	{
		//Compare only if the nucleotide is not "N"
		#ifdef NOT_COUNT_NUCLEOTIDE_N
		if(bam_seq[i] != 'N')
		#endif
		{
			//Recalibrate quality
			qual_index = bam_quals[i] - bam_info->min_qual;
			delta_r = bam_info->qual_delta[qual_index];

			matrix_index = qual_index * bam_info->num_cycles + i;
			delta_rc = bam_info->qual_cycle_delta[matrix_index];

			//dont take prev dinuc in first cycle (delta = 0)
			delta_rd = 0.0;
			if(i > 0)
			{
				recal_get_dinuc(bam_seq[i-1], bam_seq[i], &dinuc);
				if(dinuc != d_X)
				{
					matrix_index = qual_index * bam_info->num_dinuc + i;
					delta_rd = bam_info->qual_dinuc_delta[matrix_index];
				}
			}

			//Recalibration formula
			double calidad = (double)bam_quals[i];
			if(calidad >= MIN_QUALITY_TO_STAT)
			{
				double res = bam_info->total_estimated_Q + bam_info->total_delta + delta_r + delta_rc /*+ delta_rd*/;
				res_quals[i] = (char)res;
			}
			/*else
			{
				res_quals[i] = (char)calidad;
			}*/
		}
	}

	memcpy(bam1_qual(alig), res_quals, bam_seq_l);

	//Memory free
	free(bam_seq);
	free(res_quals);

	return NO_ERROR;
}

