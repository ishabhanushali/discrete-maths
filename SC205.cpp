#include<conio.h>
#include<iostream>
#include<vector>
#include<string>
#include<algorithm>
#include <windows.h>
#include <stdio.h>
using namespace std;

// Function to perform global sequence alignment using Needleman-Wunsch algorithm
void globalAlignment(const string& seq1, const string& seq2, int gapPenalty, int mismatchPenalty, int matchScore)
{
    int m = seq1.length();
    int n = seq2.length();

    vector<vector<int>> dp(m + 1, vector<int>(n + 1, 0));

    // Initialize the first row and column of the DP matrix
    for (size_t i = 0; i <= m; i++)
        dp[i][0] = i * gapPenalty;

    for (size_t j = 0; j <= n; j++)
        dp[0][j] = j * gapPenalty;

    // Fill the DP matrix
    for (size_t i = 1; i <= m; i++)
    {
        for (size_t j = 1; j <= n; j++)
        {
            int match = (seq1[i - 1] == seq2[j - 1]) ? matchScore : mismatchPenalty;
            dp[i][j] = max(dp[i - 1][j - 1] + match, max(dp[i - 1][j] + gapPenalty, dp[i][j - 1] + gapPenalty));
        }
    }

    // Traceback to find the aligned sequences
    string alignedSeq1 = "";
    string alignedSeq2 = "";
    int i = m;
    int j = n;

    while (i > 0 || j > 0)
    {
        if (i > 0 && j > 0 && dp[i][j] == dp[i - 1][j - 1] + ((seq1[i - 1] == seq2[j - 1]) ? matchScore : mismatchPenalty))
        {
            alignedSeq1 = seq1[i - 1] + alignedSeq1;
            alignedSeq2 = seq2[j - 1] + alignedSeq2;
            i--;
            j--;
        }
        else if (i > 0 && dp[i][j] == dp[i - 1][j] + gapPenalty)
        {
            alignedSeq1 = seq1[i - 1] + alignedSeq1;
            alignedSeq2 = "-" + alignedSeq2;
            i--;
        }
        else
        {
            alignedSeq1 = "-" + alignedSeq1;
            alignedSeq2 = seq2[j - 1] + alignedSeq2;
            j--;
        }
    }

    // Print the aligned sequences
    cout<< endl << "Global Alignment:" << endl;
    cout<<endl<<"Optimal Global Alignment Score: "<<dp[m][n]<<endl;
    cout << "Aligned Sequence 1: " << alignedSeq1 << endl;
    cout << "Aligned Sequence 2: " << alignedSeq2 << endl;
}

// Function to perform local sequence alignment using Smith-Waterman algorithm
void localAlignment(const string& seq1, const string& seq2, int gapPenalty, int mismatchPenalty, int matchScore)
{
    int m = seq1.length();
    int n = seq2.length();

    vector<vector<int>> dp(m + 1, vector<int>(n + 1, 0));
    int maxScore = 0;
    pair<int, int> maxPos;

    // Fill the DP matrix
    for (size_t i = 1; i <= m; i++)
    {
        for (size_t j = 1; j <= n; j++)
        {
            int match = (seq1[i - 1] == seq2[j - 1]) ? matchScore : mismatchPenalty;
            dp[i][j] = max(0, max(dp[i - 1][j - 1] + match, max(dp[i - 1][j] + gapPenalty, dp[i][j - 1] + gapPenalty)));

            // Update the maximum score and its position
            if (dp[i][j] > maxScore)
            {
                maxScore = dp[i][j];
                maxPos = make_pair(i, j);
            }
        }
    }

    // Traceback to find the aligned sequences
    string alignedSeq1 = "";
    string alignedSeq2 = "";
    int i = maxPos.first;
    int j = maxPos.second;

    while (i > 0 && j > 0 && dp[i][j] != 0)
    {
        if (dp[i][j] == dp[i - 1][j - 1] + ((seq1[i - 1] == seq2[j - 1]) ? matchScore : mismatchPenalty))
        {
            alignedSeq1 = seq1[i - 1] + alignedSeq1;
            alignedSeq2 = seq2[j - 1] + alignedSeq2;
            i--;
            j--;
        }
        else if (dp[i][j] == dp[i - 1][j] + gapPenalty)
        {
            alignedSeq1 = seq1[i - 1] + alignedSeq1;
            alignedSeq2 = "-" + alignedSeq2;
            i--;
        }
        else if (dp[i][j] == dp[i][j - 1] + gapPenalty)
        {
            alignedSeq1 = "-" + alignedSeq1;
            alignedSeq2 = seq2[j - 1] + alignedSeq2;
            j--;
        }
    }

    // Print the aligned sequences
    cout << endl <<"Local Alignment:" << endl;
    cout<<endl<<"Optimal Local Alignment Score: "<<maxScore<<endl;
    cout << "Aligned Sequence 1: " << alignedSeq1 << endl;
    cout << "Aligned Sequence 2: " << alignedSeq2 << endl;
}

int main()
{
    string seq1, seq2;

    cout <<endl<< "------------------------------GLOBAL AND LOCAL SEQUENCE ALIGNMENT------------------------------" << endl << endl << endl;
    cout << "Enter any two biological sequences to find its global and local alignment by Needleman-Wunsch Algorithm and Smith-Waterman Algorithm respectively.";
    cout << endl << endl << endl;
    cout << "Enter First Sequence: ";
    cin >> seq1;

    cout << "Enter Second Sequence: ";
    cin >> seq2;

    cout << endl <<"------------------------------------------------------------------------------------------------"<<endl;
    cout << endl << "Scoring Mechanism" << endl;

    int gapPenalty, mismatchPenalty, matchScore;

    cout << endl << "Enter Match Score: ";
    cin >> matchScore;
    cout << "Enter Mismatch Penalty: ";
    cin >> mismatchPenalty;
    cout << "Enter Gap Penalty(Indel Score): ";
    cin >> gapPenalty;

    cout << endl <<"------------------------------------------------------------------------------------------------"<<endl;
    globalAlignment(seq1, seq2, gapPenalty, mismatchPenalty, matchScore);
    localAlignment(seq1, seq2, gapPenalty, mismatchPenalty, matchScore);

    cout << endl <<"------------------------------------------------------------------------------------------------"<<endl;
    cout<<endl<<"Press Any Key to Exit.."<<endl;
    _getch();
    return 0;
}
