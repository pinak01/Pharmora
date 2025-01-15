
import {
  BellIcon,
  CalendarIcon,
  FileTextIcon,
  GlobeIcon,
  InputIcon,
} from "@radix-ui/react-icons";
import bg1 from "@/app/assets/dna.jpg"
import alz from "@/app/assets/alzdark.png"
import peptides from "@/app/assets/pepdark.png"
import { BentoCard, BentoGrid } from "@/components/ui/bento-grid";
import Image from "next/image";

const features = [
  {
    style: {
      background: `url(${alz.src})`, width: '100%',
      height: '100%', backgroundPosition: 'center', backgroundSize: 'cover'
    },
    name: "Bioactivity Predictor (pIC50)",
    description: "A tool to predict the bioactivity of a molecule.",
    href: "/bioactivity-predictor",
    cta: "Learn more",
    background: <img className="absolute -right-20 -top-20 opacity-60" />,

    className: "lg:row-start-1 lg:row-end-4 lg:col-start-2 lg:col-end-3",
  },
  {
    style: {
      background: `url(${peptides.src})`, width: '100%',
      height: '100%', backgroundPosition: 'center', backgroundSize: 'cover'
    },
    name: "Antimicrobial Activity Predictor for Peptides",
    description: "A tool to predict the antimicrobial activity of peptides.",
    href: "/",
    cta: "Learn more",
    background: <img className="absolute -right-20 -top-20 opacity-60" />,
    className: "lg:col-start-1 lg:col-end-2 lg:row-start-1 lg:row-end-3",
  },
  {
    style: {
      background: `url(${alz.src})`, width: '100%',
      height: '100%', backgroundPosition: 'center'
    },
    name: "Solubility Predictor",
    description: "A tool to predict the log(P) value of drugs",
    href: "/solubility-predictor",
    cta: "Learn more",
    background: <img className="absolute -right-20 -top-20 opacity-60" />,
    className: "lg:col-start-1 lg:col-end-2 lg:row-start-3 lg:row-end-4",
  },
  {
    style: {
      background: `url(${bg1.src})`, width: '100%',
      height: '100%', backgroundPosition: 'center', backgroundSize: 'cover'
    },
    name: "Protein Structure Predictor",
    description: "A comprehensive tool to predict protein structure.",
    href: "http://localhost:8501",
    cta: "Learn more",
    background: <img className="absolute -right-20 -top-20 opacity-60" />,
    className: "lg:col-start-3 lg:col-end-3 lg:row-start-1 lg:row-end-2",
  },
  {
    style: {
      background: `url(${alz.src})`, width: '100%',
      height: '100%', backgroundPosition: 'center', backgroundSize: 'cover'
    },
    name: "Molecular Descriptor Calculator",
    description:
      "An intuitive tool to calculate molecular descriptors.",
    href: "/",
    cta: "Learn more",
    background: <img className="absolute -right-20 -top-20 opacity-60" />,
    className: "lg:col-start-3 lg:col-end-3 lg:row-start-2 lg:row-end-4",
  },
  {
    style: {
      background: `url(${bg1.src})`, width: '100%',
      height: '100%', backgroundPosition: 'center', backgroundSize: 'cover'
    },
    name: "Lipinski’s Rule of Five Filter for Drugs",
    description:
      "A tool to filter out molecules that obey Lipinski's Rule of Five.",
    href: "/",
    cta: "Learn more",
    background: <img className="absolute -right-20 -top-20 opacity-60" />,
    className: "lg:col-start-1 lg:col-end-3 lg:row-start-4 lg:row-end-4",
  },
  {
    style: {
      background: `url(${peptides.src})`, width: '100%',
      height: '100%', backgroundPosition: 'center', backgroundSize: 'cover'
    },
    name: "New Molecule Discovery for Specific Targets",
    description:
      "A tool to discover new molecules for specific targets.",
    href: "/",
    cta: "Learn more",
    background: <img className="absolute -right-20 -top-20 opacity-60" />,
    className: "lg:col-start-3 lg:col-end-4 lg:row-start-4 lg:row-end-4",
  },
];

export async function Bento() {
  return (
    <>
      <div className="mt-24 flex items-center justify-center w-full bg-gradient-to-r from-primary-500 via-primary-600 to-primary-700">
        <h1 className="text-5xl font-bold pb-10 text-white">Explore our Tools</h1>
        </div>

      <BentoGrid className="bg-background mx-auto mt-2 w-10/12 lg:grid-rows-4">
        {features.map((feature) => (
          <BentoCard key={feature.name} {...feature} />
        ))}
      </BentoGrid>
    </>

  );
}
