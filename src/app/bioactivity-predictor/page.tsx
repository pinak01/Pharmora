import Navbar from '@/components/navbar'
import React from 'react'

const BioactivityPredictor = () => {
  return (
    <div className="relative flex size-full min-h-screen flex-col bg-[#F8F9FB] group/design-root overflow-x-hidden" style={{ fontFamily: "Inter, Noto Sans, sans-serif" }}>
      <div className="layout-container flex h-full grow flex-col">
        <Navbar></Navbar>
        <div className="px-40 flex flex-1 justify-center py-5">
          <div className="layout-content-container flex flex-col max-w-[960px] flex-1">
            <div className="flex flex-wrap justify-between gap-3 p-4">
              <p className="text-[#141C24] tracking-light text-[32px] font-bold leading-tight min-w-72">Bioactivity Predictor</p>
            </div>
            <p className="text-[#141C24] text-base font-normal leading-normal pb-3 pt-1 px-4">
              Predict the bioactivity of your molecules with our state-of-the-art model. Enter your SMILES strings below, and we will return the PubChem fingerprints and predicted
              pIC50 values.
            </p>
            <div className="flex max-w-[480px] flex-wrap items-end gap-4 px-4 py-3">
              <label className="flex flex-col min-w-40 flex-1">
                <textarea
                  placeholder="Enter SMILES strings here"
                  className="form-input flex w-full min-w-0 flex-1 resize-none overflow-hidden rounded text-[#141C24] focus:outline-0 focus:ring-0 border border-[#D4DBE8] bg-[#F8F9FB] focus:border-[#D4DBE8] min-h-36 placeholder:text-[#3F5374] p-[15px] text-base font-normal leading-normal"
                ></textarea>
              </label>
            </div>
            <div className="flex px-4 py-3">
              <button
                className="flex min-w-[84px] max-w-[480px] cursor-pointer items-center justify-center overflow-hidden rounded h-12 px-5 flex-1 bg-[#F4C753] text-[#141C24] text-base font-bold leading-normal tracking-[0.015em]"
              >
                <span className="truncate">Analyze</span>
              </button>
            </div>
            <div className="px-4 py-3 @container">
              <div className="flex overflow-hidden rounded border border-[#D4DBE8] bg-[#F8F9FB]">
                <table className="flex-1">
                  <thead>
                    <tr className="bg-[#F8F9FB]">
                      <th className="table-6de8e8e6-cbce-4671-9d40-adfa4048d8e4-column-120 px-4 py-3 text-left text-[#141C24] w-[400px] text-sm font-medium leading-normal">SMILES</th>
                      <th className="table-6de8e8e6-cbce-4671-9d40-adfa4048d8e4-column-240 px-4 py-3 text-left text-[#141C24] w-[400px] text-sm font-medium leading-normal">
                        PubChem Fingerprint
                      </th>
                      <th className="table-6de8e8e6-cbce-4671-9d40-adfa4048d8e4-column-360 px-4 py-3 text-left text-[#141C24] w-[400px] text-sm font-medium leading-normal">
                        Predicted pIC50
                      </th>
                    </tr>
                  </thead>
                  <tbody>
                    <tr className="border-t border-t-[#D4DBE8]">
                      <td className="table-6de8e8e6-cbce-4671-9d40-adfa4048d8e4-column-120 h-[72px] px-4 py-2 w-[400px] text-[#3F5374] text-sm font-normal leading-normal">
                        CC(C)(C)OC(=O)N1CCC(Oc2ccccc2)C1
                      </td>
                      <td className="table-6de8e8e6-cbce-4671-9d40-adfa4048d8e4-column-240 h-[72px] px-4 py-2 w-[400px] text-[#3F5374] text-sm font-normal leading-normal">
                        0101010000010000000000000000000000000000000000000000...
                      </td>
                      <td className="table-6de8e8e6-cbce-4671-9d40-adfa4048d8e4-column-360 h-[72px] px-4 py-2 w-[400px] text-[#3F5374] text-sm font-normal leading-normal">6.23</td>
                    </tr>
                    <tr className="border-t border-t-[#D4DBE8]">
                      <td className="table-6de8e8e6-cbce-4671-9d40-adfa4048d8e4-column-120 h-[72px] px-4 py-2 w-[400px] text-[#3F5374] text-sm font-normal leading-normal">
                        CC(C)(C)OC(=O)N1CCC(Oc2ccccc2)C1
                      </td>
                      <td className="table-6de8e8e6-cbce-4671-9d40-adfa4048d8e4-column-240 h-[72px] px-4 py-2 w-[400px] text-[#3F5374] text-sm font-normal leading-normal">
                        0101010000010000000000000000000000000000000000000000...
                      </td>
                      <td className="table-6de8e8e6-cbce-4671-9d40-adfa4048d8e4-column-360 h-[72px] px-4 py-2 w-[400px] text-[#3F5374] text-sm font-normal leading-normal">6.23</td>
                    </tr>
                  </tbody>
                </table>
              </div>
            </div>
            <div className="flex px-4 py-3 justify-end">
              <button
                className="flex min-w-[84px] max-w-[480px] cursor-pointer items-center justify-center overflow-hidden rounded h-10 px-4 bg-[#E4E9F1] text-[#141C24] text-sm font-bold leading-normal tracking-[0.015em]"
              >
                <span className="truncate">Download CSV</span>
              </button>
            </div>
            <div className="flex px-4 py-3 justify-end">
              <button
                className="flex min-w-[84px] max-w-[480px] cursor-pointer items-center justify-center overflow-hidden rounded h-10 px-4 bg-[#E4E9F1] text-[#141C24] text-sm font-bold leading-normal tracking-[0.015em]"
              >
                <span className="truncate">Download PDF</span>
              </button>
            </div>
          </div>
        </div>
      </div>
    </div>
  )
}

export default BioactivityPredictor
