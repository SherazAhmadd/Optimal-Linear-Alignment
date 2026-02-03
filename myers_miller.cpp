// The implementation of myers-miller alignment algorithm in C++
#include <iostream>  
#include <string>    
#include <vector>    
#include <algorithm> 
#include <fstream>   

using namespace std;


const int MATCH_SCORE      =  2;   
const int MISMATCH_PENALTY = -1;   
const int GAP_PENALTY      = -2;   
// ______________________________________________________________
struct AlignmentResult {
    string aligned_reference; 
    string aligned_query;     
    int total_score;          
};

string read_fasta_file(string filename);
vector<int> calculate_last_row_scores(string reference_seq, string query_seq);
AlignmentResult perform_alignment(string reference_seq, string query_seq);

string read_fasta_file(string filename) {
    ifstream input_file(filename); 
    
    if (!input_file.is_open()) {
        cout << "The fasta sequences are not correct/opened " << filename << endl;
        exit(1);
    }

    string line;
    string pure_sequence = "";

    while (getline(input_file, line)) {
     
        if (line.length() == 0) {
            continue;
        }

        if (line[0] == '>') {
            continue;
        }

        for (int i = 0; i < line.length(); i++) {
            char base = line[i];
            if (base != ' ' && base != '\n' && base != '\r') {
                pure_sequence += base;
            }
        }
    }

    input_file.close();
    return pure_sequence;
}



vector<int> calculate_last_row_scores(string reference_seq, string query_seq) {
    
    int ref_len = reference_seq.length();
    int query_len = query_seq.length();

    vector<int> previous_row(query_len + 1);
    vector<int> current_row(query_len + 1);

    for (int j = 0; j <= query_len; j++) {
        previous_row[j] = j * GAP_PENALTY;
    }

    for (int i = 1; i <= ref_len; i++) {
    
        current_row[0] = i * GAP_PENALTY;

        for (int j = 1; j <= query_len; j++) {
            
          
            char ref_base = reference_seq[i - 1];
            char query_base = query_seq[j - 1];
            
            int diagonal_score;
            if (ref_base == query_base) {
                diagonal_score = previous_row[j - 1] + MATCH_SCORE;
            } else {
                diagonal_score = previous_row[j - 1] + MISMATCH_PENALTY;
            }

            
            int up_score = previous_row[j] + GAP_PENALTY;
            int left_score = current_row[j - 1] + GAP_PENALTY;
            int max_score = diagonal_score;
            
            if (up_score > max_score) {
                max_score = up_score;
            }
            if (left_score > max_score) {
                max_score = left_score;
            }

            current_row[j] = max_score;
        }

        previous_row = current_row;
    }

    return previous_row;
}


//  The main, MYERS-MILLER recursive function

AlignmentResult perform_alignment(string reference_seq, string query_seq) {
    
    int ref_len = reference_seq.length();
    int query_len = query_seq.length();
    
    AlignmentResult result;

    if (ref_len == 0) {
        result.aligned_reference = "";
        result.aligned_query = "";
        
        for (int k = 0; k < query_len; k++) {
            result.aligned_reference += "-";
            result.aligned_query += query_seq[k];
        }
        result.total_score = query_len * GAP_PENALTY;
        return result;
    }

    if (query_len == 0) {
        result.aligned_reference = "";
        result.aligned_query = "";

        for (int k = 0; k < ref_len; k++) {
            result.aligned_reference += reference_seq[k];
            result.aligned_query += "-";
        }
        result.total_score = ref_len * GAP_PENALTY;
        return result;
    }


    if (ref_len == 1) {
        int best_score = -999999;
        int split_index = -1;

        for (int j = 0; j < query_len; j++) {
            int current_score = 0;
            
            if (reference_seq[0] == query_seq[j]) {
                current_score += MATCH_SCORE;
            } else {
                current_score += MISMATCH_PENALTY;
            }

            current_score += (query_len - 1) * GAP_PENALTY;

            if (current_score > best_score) {
                best_score = current_score;
                split_index = j;
            }
        }


        int all_gaps_score = (query_len + 1) * GAP_PENALTY;
        
        if (all_gaps_score > best_score) {
            result.aligned_reference = reference_seq + string(query_len, '-');
            result.aligned_query = "-" + query_seq;
            result.total_score = all_gaps_score;
        } else {
            result.aligned_reference = "";
            result.aligned_query = query_seq;

            for(int k=0; k < split_index; k++) result.aligned_reference += "-";
            result.aligned_reference += reference_seq[0];
            for(int k=split_index+1; k < query_len; k++) result.aligned_reference += "-";
            
            result.total_score = best_score;
        }
        return result;
    }


    int mid_point = ref_len / 2;

    string top_half_ref = reference_seq.substr(0, mid_point);
    string bottom_half_ref = reference_seq.substr(mid_point);

    vector<int> scores_top = calculate_last_row_scores(top_half_ref, query_seq);

    string bottom_half_ref_rev = bottom_half_ref;
    string query_seq_rev = query_seq;
    
    reverse(bottom_half_ref_rev.begin(), bottom_half_ref_rev.end());
    reverse(query_seq_rev.begin(), query_seq_rev.end());

    vector<int> scores_bottom = calculate_last_row_scores(bottom_half_ref_rev, query_seq_rev);

    int max_combined_score = -999999;
    int best_split_k = -1;

    for (int k = 0; k <= query_len; k++) {
        int current_total = scores_top[k] + scores_bottom[query_len - k];
        
        if (current_total > max_combined_score) {
            max_combined_score = current_total;
            best_split_k = k;
        }
    }

    string query_left_part = query_seq.substr(0, best_split_k);
    string query_right_part = query_seq.substr(best_split_k);

    AlignmentResult result_left = perform_alignment(top_half_ref, query_left_part);
    AlignmentResult result_right = perform_alignment(bottom_half_ref, query_right_part);

    result.aligned_reference = result_left.aligned_reference + result_right.aligned_reference;
    result.aligned_query = result_left.aligned_query + result_right.aligned_query;
    result.total_score = result_left.total_score + result_right.total_score;

    return result;
}


int main(int argc, char* argv[]) {
    
    if (argc != 3) {
        cout << "./myers <reference_file.fasta> <query_file.fasta>" << endl;
        return 1;
    }


    cout << "++++++++ Myers-Miller pairwise-sequence aligner ++++++++" << endl;

    string reference_file = argv[1];
    string query_file = argv[2];

    string ref_seq = read_fasta_file(reference_file);
    string query_seq = read_fasta_file(query_file);

    cout << "Reference Sequence Length: " << ref_seq.length() << " bases" << endl;
    cout << "Query Sequence Length:     " << query_seq.length() << " bases" << endl;
    cout << "Calculating optimal alignment..." << endl;
    
    AlignmentResult final_result = perform_alignment(ref_seq, query_seq);

    cout << "Alignment finished." << endl;
    cout << "Final Alignment Score: " << final_result.total_score << endl;
    
    ofstream output_file("alignment_result.txt");
    output_file << "Alignment Score: " << final_result.total_score << "\n\n";
    output_file << "Reference: " << final_result.aligned_reference << "\n";
    output_file << "Query:     " << final_result.aligned_query << "\n";
    output_file.close();

    return 0;
}