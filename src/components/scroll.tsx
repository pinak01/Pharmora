"use client";
import React from "react";
import { ContainerScroll } from "./ui/container-scroll-animation";
import Image from "next/image";
import heroImage from "@/app/assets/ml.jpg";
import solubility from "@/app/assets/solubility.jpeg";
export function HeroScrollDemo() {
  return (
    <div className="flex bg-background flex-col overflow-hidden">
      <ContainerScroll height="70" orientation="left" class2="pt-10 md:pt-40 w-full relative" className="pb-0"
        titleComponent={
          <>
            <h1 className="text-3xl font-semibold text-foreground ">
              Unleash the power of <br />
              <span className="text- w-4/6 md:text-[5rem] font-bold mt-0 leading-none">
                Artificial Intelligence
              </span>
            </h1>
          </>
        }
        textDescription="We are an all-in-one bioinformatics research platform designed to streamline drug discovery and molecular research. The platform offers advanced tools such as a Bioactivity Predictor (pIC50), New Molecule Discovery for Specific Targets, Solubility Predictor, DNA Nucleotide Counter, Antimicrobial Activity Predictor for Peptides, Molecular Descriptor Calculator, and a Lipinski’s Rule of Five Filter for Drugs. It is an ideal solution for researchers and scientists aiming to explore drug-like molecules, predict biological activity, and optimize lead compounds in a seamless, user-friendly interface."
        textHeader="Your One Stop Solution for Bioinformatics Research"
      >
        <Image
          src={heroImage}
          alt="hero"
          height={470}
          width={400}
          className="rounded-2xl object-contain h-full object-left-top"
          draggable={true}
        />
      </ContainerScroll>
      <ContainerScroll height="10" orientation="right"
      className="pt-0"
      class2="w-full relative"
        titleComponent={
          <>
          </>
        }
        textDescription="We are an all-in-one bioinformatics research platform designed to streamline drug discovery and molecular research. The platform offers advanced tools such as a Bioactivity Predictor (pIC50), New Molecule Discovery for Specific Targets, Solubility Predictor, DNA Nucleotide Counter, Antimicrobial Activity Predictor for Peptides, Molecular Descriptor Calculator, and a Lipinski’s Rule of Five Filter for Drugs. It is an ideal solution for researchers and scientists aiming to explore drug-like molecules, predict biological activity, and optimize lead compounds in a seamless, user-friendly interface."
        textHeader="Your One Stop Solution for Bioinformatics Research"
      >
        <Image
          src={solubility}
          alt="hero"
          height={470}
          width={400}
          className="rounded-2xl object-contain h-full object-left-top"
          draggable={true}
        />
      </ContainerScroll>
    </div>
  );
}
