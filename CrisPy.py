from Bio import SeqIO, pairwise2, Seq
import copy
import math
import re

# List of bases to loop over
TRACE_LIST = ['A','G','C','T']

# Simple class to parse out ab1 file format data
class ABIparse(object):
    def __init__(self, file_name):
        record = SeqIO.read(file_name,"abi")
        self.trace = {}
        self.trace['A'] = list(record.annotations['abif_raw']['DATA10'])
        self.trace['G'] = list(record.annotations['abif_raw']['DATA9'])
        self.trace['C'] = list(record.annotations['abif_raw']['DATA12'])
        self.trace['T'] = list(record.annotations['abif_raw']['DATA11'])
        self.trace['sequence'] = record.annotations['abif_raw']['PBAS2'] 
        self.trace['base_pos'] = record.annotations['abif_raw']['PLOC2']
        self.trace['seq_length'] = len(self.trace['A'])

    def get_trace(self):
        return self.trace

    def get_sequence(self):
        return self.trace['sequence'], self.trace['seq_length']

    def get_base_pos(self):
        return self.trace['base_pos']
# End ABIparse Class

# This class is a modified version of SeqDoc (Crowe, M. L. (2005). BMC Bioinformatics, 6(1), 133)
# its purpose is to normalize, align, and find the difference of two sanger sequence traces
# the script was converted to python and updated for integration with Sequalizer
class SeqDoc(object):
    def __init__(self, ref_file, test_file):
        ref_abi = ABIparse(ref_file)
        test_abi = ABIparse(test_file)      
        self.ref_trace = ref_abi.trace
        self.test_trace = test_abi.trace 

    def normalize_data(self, trace_data):

        # can only normalize for larger sequences
        if trace_data['seq_length'] < 1100:
            raise Exception('sequence too short')

        # Need to create new reference arrays, since normalization alters the existing
        # values in the hash
        orig_trace = {}
        for letter in TRACE_LIST:
            orig_trace[letter] = trace_data[letter].copy()

        # Calculate normalized datapoints, starting with the middle points
        for datapoint in range(500, trace_data['seq_length']-501):
            total_sum =0
            # Normalize to 100. Divide by sum of all values, multiply by number of  
            # values, and multiply by 100;
            # adding up the 1000 values around datapoint
            for letter in TRACE_LIST:
                temp_array = orig_trace[letter][datapoint-500:datapoint+500]
                total_sum += sum(temp_array)
            for letter in TRACE_LIST:
                # Blank sequence can cause problems through division by zero errors.
                # Deleting trailing blank sequence helps, but just put in a default value
                # for totalsum in case of problems.
                if total_sum == 0:
                    total_sum = 1000
                # Calculate normalized data
                trace_data[letter][datapoint] = int((orig_trace[letter][datapoint] / total_sum) * 4000 * 100)
        
        # Now do first 500 - special case, since can't do 500 before. Instead just 
        # take all points before. Not so critical anyway, since data quality is poor
        # at start so any mismatches will be unreliable anyway.
        # Start by initialising totalsum
        total_sum = 0
        for letter in TRACE_LIST:
            temp_array = orig_trace[letter][0:500]
            total_sum += sum(temp_array)
        # Blank sequence can cause problems through division by zero errors.
        # Deleting trailing blank sequence helps, but just put in a default value
        # for totalsum in case of problems.     
        if total_sum == 0:
            total_sum = 1000

        for datapoint in range(0,499):
            #Can do 500 after though
            end = datapoint + 500
            # Normalize to 100. Divide by sum of all values, multiply by number of
            # values, and multiply by 100;
            for letter in TRACE_LIST:
                if total_sum == 0:
                    total_sum = 1000
                # Calculate normalized data
                trace_data[letter][datapoint] = int((orig_trace[letter][datapoint] / total_sum) * end * 4 * 100)
            for letter in TRACE_LIST:
                # Add next value to totalsum, to keep 500 values after
                total_sum += orig_trace[letter][end + 1]

        # Finally the last 500 - again a special case, as can't do 500 after. Instead 
        # just take all points after. Not so critical anyway, since data quality is
        # poor at end so any mismatches will be unreliable anyway.
        # Again start by initialising totalsum
        total_sum = 0
        for letter in TRACE_LIST:
            temp_array = orig_trace[letter][trace_data['seq_length'] - 1000: trace_data['seq_length'] - 1]
            total_sum += sum(temp_array)
        # Identify start and finish points
        (first, last) = (trace_data['seq_length']-500, trace_data['seq_length'] - 1)
        for datapoint in range(first, last):
            start = datapoint - 500
            # Normalize to 100. Divide by sum of all values, multiply by number of
            # values, and multiply by 100;
            for letter in TRACE_LIST:
                if total_sum == 0:
                    total_sum = 1000
                # Calculate normalized data
                trace_data[letter][datapoint] = int((orig_trace[letter][datapoint] / total_sum) * (last-start) * 4 *100)
            # Subtract first value from totalsum, to keep to 500 point before test
            for letter in TRACE_LIST:
                total_sum -= abs(orig_trace[letter][start])

    def get_best_align(self, ref_trace, test_trace):
        # This does an alignment of the first 1000 datapoints using a range of offsets
        # from -200 to 200. The best alignment is picked on the basis of having the
        # lowest score from datapoint 200 to 1000, and is used to allow for any
        # variation in start position of the two sequences.
        scores = {}
        temp_ref = {}
        temp_test = {}
        for offset in range(-200, 200, 20):
            # Create temporary hashes, since a hash reference is passed to the function
            # and otherwise the real hash will be modified
            for letter in TRACE_LIST:
                temp_ref[letter] = ref_trace[letter].copy()
                temp_test[letter] = test_trace[letter].copy()
            # Do a partial alignment (first 1000 datapoints)
            temp_ref, temp_test = self._align(temp_ref, temp_test, offset, 1000)
            # Work out the score for that alignment
            scores[offset] = self._get_score(200, 1000, 0, temp_ref, temp_test)
        # Sort the scores to find out the lowest, and record the value of that offset
        offset = sorted(scores.items(), key=lambda x:x[1])
        # Once the best alignment has been determined, then do it for real
        self._align(ref_trace, test_trace, offset[0][0], len(test_trace['A'])+offset[0][0])

    def _align(self, ref, test, min_index, trace_length):
        # This takes the normalized traces and returns a best alignment of the two.
        # Rows are added to or deleted from the test trace to keep the alignment.
        # Best alignment is calculated by minimising the difference between the test
        # and reference sequence over the next 30 datapoints. It is adjusted every five
        # bases. Inserted lines are just a duplicate of the previous line.

        # Add/delete the appropriate number of lines to the test sequence to correct
        # for the offset value
        if min_index < 0:
            # Need to add lines to test sequence, since it's behind
            for i in range(min_index, 0):
                for letter in TRACE_LIST:
                    test[letter].insert(0,test[letter][0])
        elif min_index > 0:
            # Need to delete lines from test sequence, since it's ahead
            for letter in TRACE_LIST:
                del test[letter][:min_index]
        # Make a note of the offset value for datapoint numbering
        test['initial_offset'] = min_index
        ref['initial_offset'] = 0

        # Now check alignments
        for i in range(0, trace_length-1):
            # Each third entry (starting from 1), check alignment
            if i%3:
                continue
            start_pos = i
            end_pos = i + 30
            # Compare the scores in the current alignment with those one data point in 
            # either direction
            score = self._get_score(start_pos, end_pos, 0, ref, test)
            pre_score = self._get_score(start_pos, end_pos, -1, ref, test)
            post_score = self._get_score(start_pos, end_pos, 1, ref, test)
            if (score == 'no_score') or (pre_score == 'no_score') or (post_score == 'no_score'):
                break
            # Work out offset
            # Default is 0; score is the lowest of the three
            if (score < pre_score) and (score < post_score):
                offset = 0
            # if pre-score is the lowest, then offset = -1
            elif (pre_score < score) and (pre_score < post_score):
                offset = -1
            # if post-score is the lowest, then offset = 1
            elif (post_score < pre_score) and (post_score < score):
                offset = 1
            # If in doubt, default to no change
            else:
                offset = 0
            # Now insert or delete lines as required
            if offset == 1:
                # The reference sample is behind, need to delete a row from test
                for letter in TRACE_LIST:
                    del test[letter][i]
            elif offset == -1:
                # The reference sample is ahead, need to add a row to test
                for letter in TRACE_LIST:
                    test[letter].insert(i, test[letter][i])

        # Reset length of test sequence to correct value
        test['seq_length'] = len(test['A'])

        return ref, test

    def _get_score(self, start, end, offset, ref, test):
        # Subroutine used in alignment testing - it gets the total difference between
        # the two submitted sections of array and returns it.
        score = 0
        for i in range(start, end):
            try:
                for letter in TRACE_LIST:
                    score += abs(ref[letter][i] - test[letter][i + offset])
            except IndexError:
                return 'no_score'
        return score
    
    def differences(self, ref, test):
        # Takes the two traces and calculates the difference between the two. Then 
        # squares this difference to highlight the larger (relevent) changes.
        # Also returns the length of the shortest sequence, so both are aligned on 
        # image generation
        min_index = min(test['seq_length'], ref['seq_length'])

        # Create a hash for the difference values
        diffs = {}
        for letter in TRACE_LIST:
            diffs[letter] = []
        # Loop through the traces and calculate the differences
        for i in range(0, min_index):
            # Need to do all four traces
            for letter in TRACE_LIST:
                # Get the difference
                diff = ref[letter][i] - test[letter][i]
                # Is the difference positive or negative (need to record, otherwise will 
                # be lost on squaring)
                sign_func = lambda a: (a>0) - (a<0)
                sign = sign_func(diff)
                # Stop saturation by cutting off to a max value
                if abs(diff) > 5000:
                    diff = 5000 * sign
                # Highlight differences by multiplying value by total values of OTHER
                # channels
                diffs[letter].append(diff)

            # Have now got difference for all four traces. Can accentuate real diffs
            for letter in TRACE_LIST:
                diff = diffs[letter][i]
                # Is the difference positive or negative (need to record, otherwise will 
                # be lost on squaring)
                sign_func = lambda a: (a>0) - (a<0)
                sign = sign_func(diff)
                otherchannels = 1
                # Sum all values in the other channels which have the opposite sign
                for channel in TRACE_LIST:
                    if channel == letter:
                        continue
                    value = diffs[channel][i]
                    # Ignore if the sign is the same
                    if value * sign > 0:
                        continue
                    otherchannels += value
                finaldiff = (sign * diff * diff * math.sqrt(abs(otherchannels))) / 5000
                if abs(finaldiff) > 5000:
                    finaldiff = sign * 5000
                diffs[letter][i] = finaldiff

        return min_index, diffs

    def get_all_data(self):
        # normalize data
        self.normalize_data(self.ref_trace)
        self.normalize_data(self.test_trace)
        # align the two sequences
        self.get_best_align(self.ref_trace, self.test_trace)
        # get differece traces
        align_length, diffs = self.differences(self.ref_trace, self.test_trace)
        return align_length, diffs
# End SeqDoc class

# This class uses Timothy K. Lu lab's sequalizer formula to caclulate point mutation frequency
# converted to python and modified for integration with SeqDoc and for off-target analysis
class Sequalizer(object):
    def __init__(self, ref_data, test_data, diff_data, target_sequence, target_range):
        self.ref_data = ref_data
        self.test_data = test_data
        self.diff_data = diff_data
        self.target_sequence = target_sequence
        self.target_range = target_range
        self.align_start = None
        self.align_end = None
        self.mutation_freq = []

    def _align_seqs(self, match_override):
        if match_override is not None:
            # uses an off-target alignment, on either sense or antisense strand
            if 0 < match_override < len(self.diff_data['A']):
                # sense strand
                self.align_start = match_override - 20
                self.align_end = match_override
            elif match_override < -(len(self.ref_data['A']) - len(self.diff_data['A'])):
                # antisense strand
                self.align_start = (len(self.ref_data['A']) + match_override) + 3
                self.align_end = self.align_start + 20
        else:
            # uses pairwise2 module to find best target<->sequence alignment, records start and stop index
            alignments = pairwise2.align.localms(self.ref_data['sequence'], self.target_sequence, 2, -1, -1, -0.1)
            self.align_start = alignments[0][3]
            self.align_end = alignments[0][4]

    def get_mutation_freq(self, match_override=None):
        # finds where the target sequence occurs in the reference trace
        self.mutation_freq = []
        self._align_seqs(match_override)

        # loops through each base index specified in the target sequence
        for base_index in self.target_range:
            # specifies the base type (A,G,C,T) and sequence position from target base-index
            if match_override is not None:
                if match_override > 0:
                    current_position = self.align_start + base_index
                    target_base = self.ref_data['sequence'][current_position]
                else:
                    current_position = self.align_end - base_index
                    target_base = self.ref_data['sequence'][current_position]
            else:
                current_position = self.align_start + base_index
                target_base = self.target_sequence[base_index]
            # calculates the muation frequency for each target base
            if target_base == 'G':
                trace_range = [self.ref_data['base_pos'][current_position]-2, self.ref_data['base_pos'][current_position]+2]
                g_diff = max(self.diff_data['G'][trace_range[0]:trace_range[1]])
                a_diff = min(self.diff_data['A'][trace_range[0]:trace_range[1]])
                g_ref = max(self.ref_data['G'][trace_range[0]:trace_range[1]])
                self.mutation_freq.append(round(math.sqrt(abs((g_diff-a_diff)/(8*g_ref))),3))
            elif target_base == 'C':
                trace_range = [self.ref_data['base_pos'][current_position]-2, self.ref_data['base_pos'][current_position]+2]
                c_diff = max(self.diff_data['C'][trace_range[0]:trace_range[1]])
                t_diff = min(self.diff_data['T'][trace_range[0]:trace_range[1]])
                c_ref = max(self.ref_data['C'][trace_range[0]:trace_range[1]])
                self.mutation_freq.append(round(math.sqrt(abs((c_diff-t_diff)/(8*c_ref))),3))
            else:
                self.mutation_freq.append(0)
        return self.mutation_freq
# End Sequalizer Class

# This class finds and ranks all possible target sites in a sanger sequence
# ranking system is based off of the model from the Howard M. Salis Lab
# parameters come from the UBC 2017 iGEM team
class OfftargetFinder(object):
    def __init__(self, ref_data, target_sequence):
        self.pos_weights = (0.554111551727719,0.999999999999958,0.999859588152223,0.999997460325925,0.414113900546951,
                            0.999495056671895,0.0220208959410121,0.589953049071977,0.324385855364402,2.26201959820539e-06,
                            0.0825699665148698,0.0890566149734565,0.234751499655325,3.77820298630600e-14,0.214631126793305,
                            7.04574142003494e-06,0.156869096216009,0.129156230982504,0.0428145130615625,1.58135744395507e-05,
                            0.100000000000000)
        self.ref_data = ref_data
        sequence = Seq.Seq(''.join(ref_data['sequence']))
        self.rev_ref_sequence = str(sequence.reverse_complement())
        self.target_sequence = target_sequence        

    def _match_ngg(self):
        # finds start indexes of each PAM site match
        match_indexes = []
        # match all PAM's on sense strand
        matches = re.finditer(r"\wGG",self.ref_data['sequence'])
        for match in matches:
            match_indexes.append(match.span(0)[0])
        # match all PAM's on antisense strand (labels start index as negative)
        matches = re.finditer(r"\wGG",self.rev_ref_sequence)
        for match in matches:
            match_indexes.append(-match.span(0)[0])
        return match_indexes

    def _calc_binding(self, match_indexes):
        # scores possibles sites by approx. gibson energies
        match_dict = {}
        norm_factor = sum(self.pos_weights)
        for match_index in match_indexes:
            score = 0
            start = match_index
            if start > 0:
                for i in range(0, len(self.target_sequence)):
                    position = len(self.target_sequence) - i
                    if self.target_sequence[i] == self.ref_data['sequence'][start - position]:
                        score += 0
                    else:
                        score += self.pos_weights[position]
                score = score/norm_factor
                match_dict[score] = start
            elif start < 0:
                for i in range(0, len(self.target_sequence)):
                    position = len(self.target_sequence) - i
                    if self.target_sequence[i] == self.rev_ref_sequence[start - position]:
                        score += 0
                    else:
                        score += self.pos_weights[position]
                score = score/norm_factor
                match_dict[score] = start
        return match_dict

    def get_targets(self):
        match_indexes = self._match_ngg()
        match_dict = self._calc_binding(match_indexes)
        for k,v in match_dict.copy().items():
            print("score is: ")
            print(k)
            print("starting at:")
            print(v)
            if (k > .6) or (k == 0):
                print("deleted")
                del match_dict[k]
            print("\n")
        return match_dict
# End OfftargetFinder Class
