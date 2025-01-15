'use client'

import React from 'react'

import { usePathname } from 'next/navigation'
const Navbar = () => {
    let route:string = usePathname().slice(1).toLocaleLowerCase()
    console.log(route=='solubility-predictor');
    
    return (
        <nav className="bg-[#333333] shadow">
            <div className="container flex items-center justify-end w-12/12 max-w-full p-4 mx-auto text-gray-300 capitalize">
                <a href="/bioactivity-predictor" className={`border-b-2 mx-1.5 sm:mx-6 ${route=='bioactivity-predictor' ? 'border-[#F4C753] text-gray-200' : 'hover:border-[#F4C753] border-transparent hover:text-gray-200'}`}>Bioactivity Predictor</a>
                <a href="/dna-nucleotide-counter" className={`border-b-2 hover:text-gray-200 mx-1.5 sm:mx-6 ${route=='dna-nucleotide-counter' ? 'border-[#F4C753] text-gray-200' : 'hover:border-[#F4C753] border-transparent hover:text-gray-200'}`}>DNA Nucleotide Counter</a>
                <a href="/lipinski" className={`border-b-2 hover:text-gray-200 mx-1.5 sm:mx-6 ${route=='lipinski' ? 'border-[#F4C753] text-gray-200' : 'border-transparent hover:border-[#F4C753] hover:text-gray-200'}`}>Lipinski's Rule</a>
                <a href="/molecular-descriptor" className={`border-b-2 hover:text-gray-200 mx-1.5 sm:mx-6 ${route=='molecular-descriptor' ? 'border-[#F4C753] text-gray-200' : 'border-transparent hover:border-[#F4C753] hover:text-gray-200'}`}>Molecular Descriptor</a>
                <a href="/molecule-discovery" className={`border-b-2 hover:text-gray-200 mx-1.5 sm:mx-6 ${route=='molecular-discovery' ? 'border-[#F4C753] text-gray-200' : 'border-transparent hover:border-[#F4C753] hover:text-gray-200'}`}>New Drug</a>
                <a href="/protein-structure" className={`border-b-2 hover:text-gray-200 mx-1.5 sm:mx-6 ${route=='protein-structure' ? 'border-[#F4C753] text-gray-200' : 'border-transparent hover:border-[#F4C753] hover:text-gray-200'}`}>Protein Structure Generator</a>
                <a href="/antimicrobial-detection" className={`border-b-2 hover:text-gray-200 mx-1.5 sm:mx-6 ${route=='antimicrobial-detection' ? 'border-[#F4C753] text-gray-200' : 'border-transparent hover:border-[#F4C753] hover:text-gray-200'}`}>Peptide Activity</a>
                <a href="/solubility-predictor" className={`border-b-2 hover:text-gray-200 mx-1.5 sm:mx-6 ${route=='solubility-predictor' ? 'border-[#F4C753] text-gray-200' : 'border-transparent hover:border-[#F4C753] hover:text-gray-200'}`}>Solubility Predictor</a>
            </div>
        </nav>
    )
}

export default Navbar