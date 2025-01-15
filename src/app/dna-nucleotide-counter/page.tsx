'use client'
import Navbar from '@/components/navbar'
import React, { useState } from 'react' // Import useState

// Define the state type
type OutputState = {
    counts: { A: number; T: number; G: number; C: number; N: number; };
    totalBases: number;
    gcContent: string;
} | null;

const DnaNuc = () => {
    const [output, setOutput] = useState<OutputState>(null); // State to hold output

    const handleCountBases = () => {
        // Sample DNA sequence for demonstration
        const sampleDna = "ATGCGATACGCTAGCTAGCTAGCTAGC";
        const counts = {
            A: (sampleDna.match(/A/g) || []).length,
            T: (sampleDna.match(/T/g) || []).length,
            G: (sampleDna.match(/G/g) || []).length,
            C: (sampleDna.match(/C/g) || []).length,
            N: 0 // Assuming no ambiguous bases for simplicity
        };
        const totalBases = counts.A + counts.T + counts.G + counts.C;
        const gcContent = ((counts.G + counts.C) / totalBases) * 100;

        // Set output with detailed results
        setOutput((prevState: OutputState) => ({
            ...(prevState || {}), // Use an empty object if prevState is null
            counts,
            totalBases,
            gcContent: gcContent.toFixed(2) + '%'
        }));
    };

    return (
        <div><div className="relative flex size-full min-h-screen flex-col bg-[#1A1A1A] group/design-root overflow-x-hidden">
            <div className="layout-container flex h-full grow flex-col">
                <Navbar></Navbar>
                <div className="px-40 flex flex-1 justify-center py-5">
                    <div className="layout-content-container flex flex-col w-[512px] max-w-[512px] py-5 max-w-[960px] flex-1">
                        <h2 className="text-[#FFFFFF] tracking-light text-[28px] font-bold leading-tight px-4 text-center pb-3 pt-5">Nucleotide Base Counter</h2>
                        <p className="text-[#E0E0E0] text-base font-normal leading-normal pb-3 pt-1 px-4 text-center">
                            Input a DNA sequence and get a count of nucleotide bases with a bar chart representation, GC content percentage, batch processing of FASTA/TXT files, and options to
                            download results as CSV or PDF.
                        </p>
                        <div className="flex max-w-[480px] flex-wrap items-end gap-4 px-4 py-3">
                            <label className="flex flex-col min-w-40 flex-1">
                                <textarea
                                    placeholder="Paste or type DNA sequence here"
                                    className="form-input flex w-full min-w-0 flex-1 resize-none overflow-hidden rounded text-[#FFFFFF] focus:outline-0 focus:ring-0 border border-[#444444] bg-[#2A2A2A] focus:border-[#FFD700] hover:border-[#FFD700] min-h-36 placeholder:text-[#B0B0B0] p-[15px] text-base font-normal leading-normal" // Added hover and focus border styles
                                ></textarea>
                            </label>
                        </div>
                        <div className="flex justify-stretch">
                            <div className="flex flex-1 gap-3 flex-wrap px-4 py-3 justify-end">
                                <button
                                    className="flex min-w-[84px] max-w-[480px] cursor-pointer items-center justify-center overflow-hidden rounded h-10 px-4 bg-[#444444] text-[#FFFFFF] text-sm font-bold leading-normal tracking-[0.015em]"
                                >
                                    <span className="truncate">Upload File</span>
                                </button>
                                <button
                                    className="flex min-w-[84px] max-w-[480px] cursor-pointer items-center justify-center overflow-hidden rounded h-10 px-4 bg-[#F4C753] text-[#141C24] text-sm font-bold leading-normal tracking-[0.015em]"
                                    onClick={handleCountBases} // Attach click handler
                                >
                                    <span className="truncate">Count Bases</span>
                                </button>
                            </div>
                        </div>
                        <h3 className="text-[#FFFFFF] text-lg font-bold leading-tight tracking-[-0.015em] px-4 pb-2 pt-4">Results</h3>
                        <div className="flex flex-wrap gap-4 px-4 py-6">
                            <div className="flex min-w-72 flex-1 flex-col gap-2 rounded border border-[#444444] p-6">
                                <p className="text-[#E0E0E0] text-base font-medium leading-normal">Count of Nucleotide Bases</p>
                                <div className="grid min-h-[180px] gap-x-4 gap-y-6 grid-cols-[auto_1fr] items-center py-3">
                                    <p className="text-[#B0B0B0] text-[13px] font-bold leading-normal tracking-[0.015em]">A</p>
                                    <div className="h-full flex-1"><div className="border-[#B0B0B0] bg-[#3A3A3A] border-r-2 h-full" style={{width: "90%"}}></div></div>
                                    <p className="text-[#B0B0B0] text-[13px] font-bold leading-normal tracking-[0.015em]">T</p>
                                    <div className="h-full flex-1"><div className="border-[#B0B0B0] bg-[#3A3A3A] border-r-2 h-full" style={{width: "60%"}}></div></div>
                                    <p className="text-[#B0B0B0] text-[13px] font-bold leading-normal tracking-[0.015em]">G</p>
                                    <div className="h-full flex-1"><div className="border-[#B0B0B0] bg-[#3A3A3A] border-r-2 h-full" style={{width: "30%"}}></div></div>
                                    <p className="text-[#B0B0B0] text-[13px] font-bold leading-normal tracking-[0.015em]">C</p>
                                    <div className="h-full flex-1"><div className="border-[#B0B0B0] bg-[#3A3A3A] border-r-2 h-full" style={{width: "60%"}}></div></div>
                                    <p className="text-[#B0B0B0] text-[13px] font-bold leading-normal tracking-[0.015em]">N</p>
                                    <div className="h-full flex-1"><div className="border-[#B0B0B0] bg-[#3A3A3A] border-r-2 h-full" style={{width: "100%"}}></div></div>
                                </div>
                            </div>
                        </div>
                        <div className="p-4 grid grid-cols-2">
                            <div className="flex flex-col gap-1 border-t border-solid border-t-[#444444] py-4 pr-2">
                                <p className="text-[#B0B0B0] text-sm font-normal leading-normal">GC Content</p>
                                <p className="text-[#E0E0E0] text-sm font-normal leading-normal">50.00%</p>
                            </div>
                            <div className="flex flex-col gap-1 border-t border-solid border-t-[#444444] py-4 pl-2">
                                <p className="text-[#B0B0B0] text-sm font-normal leading-normal">Total Bases</p>
                                <p className="text-[#E0E0E0] text-sm font-normal leading-normal">100</p>
                            </div>
                            <div className="flex flex-col gap-1 border-t border-solid border-t-[#444444] py-4 pr-2">
                                <p className="text-[#B0B0B0] text-sm font-normal leading-normal">A</p>
                                <p className="text-[#E0E0E0] text-sm font-normal leading-normal">25</p>
                            </div>
                            <div className="flex flex-col gap-1 border-t border-solid border-t-[#444444] py-4 pl-2">
                                <p className="text-[#B0B0B0] text-sm font-normal leading-normal">T</p>
                                <p className="text-[#E0E0E0] text-sm font-normal leading-normal">25</p>
                            </div>
                            <div className="flex flex-col gap-1 border-t border-solid border-t-[#444444] py-4 pr-2">
                                <p className="text-[#B0B0B0] text-sm font-normal leading-normal">G</p>
                                <p className="text-[#E0E0E0] text-sm font-normal leading-normal">25</p>
                            </div>
                            <div className="flex flex-col gap-1 border-t border-solid border-t-[#444444] py-4 pl-2">
                                <p className="text-[#B0B0B0] text-sm font-normal leading-normal">C</p>
                                <p className="text-[#E0E0E0] text-sm font-normal leading-normal">25</p>
                            </div>
                            <div className="flex flex-col gap-1 border-t border-solid border-t-[#444444] py-4 pr-2 col-span-2 pr-[50%]">
                                <p className="text-[#B0B0B0] text-sm font-normal leading-normal">N</p>
                                <p className="text-[#E0E0E0] text-sm font-normal leading-normal">0</p>
                            </div>
                        </div>
                        <div className="flex justify-stretch">
                            <div className="flex flex-1 gap-3 flex-wrap px-4 py-3 justify-end">
                                <button
                                    className="flex min-w-[84px] max-w-[480px] cursor-pointer items-center justify-center overflow-hidden rounded h-10 px-4 bg-[#444444] text-[#FFFFFF] text-sm font-bold leading-normal tracking-[0.015em]"
                                >
                                    <span className="truncate">Download as CSV</span>
                                </button>
                                <button
                                    className="flex min-w-[84px] max-w-[480px] cursor-pointer items-center justify-center overflow-hidden rounded h-10 px-4 bg-[#F4C753] text-[#141C24] text-sm font-bold leading-normal tracking-[0.015em]"
                                >
                                    <span className="truncate">Download as PDF</span>
                                </button>
                            </div>
                        </div>
                        {/* Display output if available */}
                        {output && (
                            <div className="output-container p-4">
                                <h3 className="text-[#FFFFFF] text-lg font-bold">Detailed Output</h3>
                                <p className="text-[#E0E0E0]">Total Bases: {output.totalBases}</p>
                                <p className="text-[#E0E0E0]">GC Content: {output.gcContent}</p>
                                <p className="text-[#E0E0E0]">Counts:</p>
                                <ul className="text-[#E0E0E0]">
                                    <li>A: {output.counts.A}</li>
                                    <li>T: {output.counts.T}</li>
                                    <li>G: {output.counts.G}</li>
                                    <li>C: {output.counts.C}</li>
                                    <li>N: {output.counts.N}</li>
                                </ul>
                            </div>
                        )}
                    </div>
                </div>
            </div>
        </div>
        </div>
    )
}

export default DnaNuc