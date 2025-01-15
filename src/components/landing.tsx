import React from "react";

import Globe from "@/components/ui/globe";
import localFont from "next/font/local";

const akeila = localFont({ src: "../fonts/Akeila.otf" });

const Landing = () => {
  return (
    <div className="z-50 w-full h-full relative flex size-full justify-center overflow-hidden bg-background md:shadow-xl">
      <div className="mt-16 z-50 block text-center justify-center w-9/12 gap-0 grid-rows-6">
        <p
          className={
            akeila.className +
            ` z-50 h-32 inline-block whitespace-pre-line text-transparent pointer-events-none bg-gradient-to-b from-neutral-600 to-white bg-clip-text text-center text-9xl font-semibold leading-none `
          }
        >
          Pharmora
        </p>
        <p className="text-transparent inline-block pointer-events-none bg-gradient-to-b from-neutral-600 to-white bg-clip-text z-50 h-10 text-center font-semibold leading-none whitespace-pre-line text-[24px] mx-auto">
          Accelerating Innovation in Bioinformatics
        </p>
      </div>
      <Globe className="top-44 z-30" />
      <div className="pointer-events-none absolute inset-0 h-full bg-[radial-gradient(circle_at_50%_200%,rgba(0,0,0,0.2),rgba(255,255,255,0))]" />
    </div>
  );
};

export default Landing;
