import Navbar from '@/components/navbar'
import React from 'react'

const SolubilityPredictor = () => {
    const css = `
    @container(max-width:120px){.table-17198cd9-e0b1-42c5-8b0d-3bbb3865882c-column-120{display: none;}}
                @container(max-width:240px){.table-17198cd9-e0b1-42c5-8b0d-3bbb3865882c-column-240{display: none;}}
                @container(max-width:360px){.table-17198cd9-e0b1-42c5-8b0d-3bbb3865882c-column-360{display: none;}}`
    return (
        <div className="relative flex size-full min-h-screen flex-col bg-[#1E1E1E] group/design-root overflow-x-hidden" style={{ fontFamily: 'Inter, "Noto Sans", sans-serif', color: '#FFFFFF' }}>
            <div className="layout-container flex h-full grow flex-col">
                <Navbar></Navbar>
                <div className="px-40 flex flex-1 justify-center py-5">
                    <div className="layout-content-container flex flex-col max-w-[960px] flex-1">
                        <div className="flex flex-wrap justify-between gap-3 p-4">
                            <div className="flex min-w-72 flex-col gap-3">
                                <p className="text-[#FFFFFF] tracking-light text-[32px] font-bold leading-tight">Solubility Predictor</p>
                                <p className="text-[#B0B0B0] text-sm font-normal leading-normal">Predict the solubility (logS) of a molecule.</p>
                            </div>
                        </div>
                        <div className="flex max-w-[480px] flex-wrap items-end gap-4 px-4 py-3">
                            <label className="flex flex-col min-w-40 flex-1">
                                <textarea
                                    placeholder="Enter SMILES (or InChI) here, one per line. Max 50 lines."
                                    className="form-input flex w-full min-w-0 flex-1 resize-none overflow-hidden rounded text-[#FFFFFF] focus:outline-0 focus:ring-0 border border-[#444444] bg-[#2A2A2A] focus:border-[#666666] min-h-36 placeholder:text-[#B0B0B0] p-[15px] text-base font-normal leading-normal hover:border-[#F4C753] transition-colors duration-300"
                                ></textarea>
                            </label>
                        </div>
                        <div className="flex px-4 py-3 justify-end">
                            <button
                                className="flex min-w-[84px] max-w-[480px] cursor-pointer items-center justify-center overflow-hidden rounded h-10 px-4 bg-[#F4C753] text-[#141C24] text-sm font-bold leading-normal tracking-[0.015em] hover:bg-[#FFD700] transition-colors duration-300"
                            >
                                <span className="truncate">Predict Solubility</span>
                            </button>
                        </div>
                        <h2 className="text-[#FFFFFF] text-[22px] font-bold leading-tight tracking-[-0.015em] px-4 pb-3 pt-5">Results</h2>
                        <div className="px-4 py-3 @container">
                            <div className="flex overflow-hidden rounded border border-[#444444] bg-[#2A2A2A] hover:shadow-lg transition-shadow duration-300">
                                <table className="flex-1">
                                    <thead>
                                        <tr className="bg-[#2A2A2A]">
                                            <th className="table-17198cd9-e0b1-42c5-8b0d-3bbb3865882c-column-120 px-4 py-3 text-left text-[#FFFFFF] w-[400px] text-sm font-medium leading-normal">SMILES</th>
                                            <th className="table-17198cd9-e0b1-42c5-8b0d-3bbb3865882c-column-240 px-4 py-3 text-left text-[#FFFFFF] w-[400px] text-sm font-medium leading-normal">logS</th>
                                            <th className="table-17198cd9-e0b1-42c5-8b0d-3bbb3865882c-column-360 px-4 py-3 text-left text-[#FFFFFF] w-60 text-sm font-medium leading-normal">Error</th>
                                        </tr>
                                    </thead>
                                    <tbody>
                                        <tr className="border-t border-t-[#444444] hover:bg-[#3A3A3A] transition-colors duration-300">
                                            <td className="table-17198cd9-e0b1-42c5-8b0d-3bbb3865882c-column-120 h-[72px] px-4 py-2 w-[400px] text-[#B0B0B0] text-sm font-normal leading-normal">
                                                C1=CC=CC=C1
                                            </td>
                                            <td className="table-17198cd9-e0b1-42c5-8b0d-3bbb3865882c-column-240 h-[72px] px-4 py-2 w-[400px] text-[#B0B0B0] text-sm font-normal leading-normal">-1.2</td>
                                            <td className="table-17198cd9-e0b1-42c5-8b0d-3bbb3865882c-column-360 h-[72px] px-4 py-2 w-60 text-sm font-normal leading-normal">
                                                <button
                                                    className="flex min-w-[84px] max-w-[480px] cursor-pointer items-center justify-center overflow-hidden rounded h-8 px-4 bg-[#444444] text-[#FFFFFF] text-sm font-medium leading-normal w-full hover:bg-[#555555] transition-colors duration-300"
                                                >
                                                    <span className="truncate">0.1</span>
                                                </button>
                                            </td>
                                        </tr>
                                        <tr className="border-t border-t-[#444444] hover:bg-[#3A3A3A] transition-colors duration-300">
                                            <td className="table-17198cd9-e0b1-42c5-8b0d-3bbb3865882c-column-120 h-[72px] px-4 py-2 w-[400px] text-[#B0B0B0] text-sm font-normal leading-normal">
                                                C1=CC=CC=C1
                                            </td>
                                            <td className="table-17198cd9-e0b1-42c5-8b0d-3bbb3865882c-column-240 h-[72px] px-4 py-2 w-[400px] text-[#B0B0B0] text-sm font-normal leading-normal">-1.2</td>
                                            <td className="table-17198cd9-e0b1-42c5-8b0d-3bbb3865882c-column-360 h-[72px] px-4 py-2 w-60 text-sm font-normal leading-normal">
                                                <button
                                                    className="flex min-w-[84px] max-w-[480px] cursor-pointer items-center justify-center overflow-hidden rounded h-8 px-4 bg-[#444444] text-[#FFFFFF] text-sm font-medium leading-normal w-full hover:bg-[#555555] transition-colors duration-300"
                                                >
                                                    <span className="truncate">0.1</span>
                                                </button>
                                            </td>
                                        </tr>
                                    </tbody>
                                </table>
                            </div>
                            <style>
                                {css}
                            </style>
                        </div>
                        <h2 className="text-[#FFFFFF] text-[22px] font-bold leading-tight tracking-[-0.015em] px-4 pb-3 pt-5">Batch Processing</h2>
                        <p className="text-[#B0B0B0] text-base font-normal leading-normal pb-3 pt-1 px-4">
                            Upload a CSV file with a column of SMILES (or InChI) to predict solubility for each molecule.
                        </p>
                        <div className="flex px-4 py-3 justify-start">
                            <button
                                className="flex min-w-[84px] max-w-[480px] cursor-pointer items-center justify-center overflow-hidden rounded h-10 px-4 bg-[#444444] text-[#FFFFFF] text-sm font-bold leading-normal tracking-[0.015em] hover:bg-[#555555] transition-colors duration-300"
                            >
                                <span className="truncate">Upload CSV</span>
                            </button>
                        </div>
                        <p className="text-[#B0B0B0] text-base font-normal leading-normal pb-3 pt-1 px-4">Download results as a CSV or PDF file.</p>
                        <div className="flex px-4 py-3 justify-start">
                            <button
                                className="flex min-w-[84px] max-w-[480px] cursor-pointer items-center justify-center overflow-hidden rounded h-10 px-4 bg-[#444444] text-[#FFFFFF] text-sm font-bold leading-normal tracking-[0.015em] hover:bg-[#555555] transition-colors duration-300"
                            >
                                <span className="truncate">Download CSV</span>
                            </button>
                        </div>
                        <div className="flex px-4 py-3 justify-start">
                            <button
                                className="flex min-w-[84px] max-w-[480px] cursor-pointer items-center justify-center overflow-hidden rounded h-10 px-4 bg-[#444444] text-[#FFFFFF] text-sm font-bold leading-normal tracking-[0.015em] hover:bg-[#555555] transition-colors duration-300"
                            >
                                <span className="truncate">Download PDF</span>
                            </button>
                        </div>
                        <h2 className="text-[#FFFFFF] text-[22px] font-bold leading-tight tracking-[-0.015em] px-4 pb-3 pt-5">Interactive Solubility Chart</h2>
                        <p className="text-[#B0B0B0] text-base font-normal leading-normal pb-3 pt-1 px-4">
                            The chart below shows the solubility of common functional groups. The solubility value is logS, where a positive value means more soluble and a negative value means
                            less soluble.
                        </p>
                    </div>
                </div>
            </div>
        </div>
    )
}

export default SolubilityPredictor