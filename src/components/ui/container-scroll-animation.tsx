"use client";
import React, { useRef } from "react";
import { useScroll, useTransform, motion, MotionValue } from "framer-motion";
import { div } from "framer-motion/client";

export const ContainerScroll = ({
  orientation,
  titleComponent,
  children,
  textDescription,
  textHeader,
  className,
  class2,
  height,
}: {
  orientation: "left" | "right";
  titleComponent: string | React.ReactNode;
  children: React.ReactNode;
  textDescription: string;
  textHeader: string;
  className?: string;
  class2?: string;
  height: string;
}) => {
  const containerRef = useRef<any>(null);
  const { scrollYProgress } = useScroll({
    target: containerRef,
  });
  const [isMobile, setIsMobile] = React.useState(false);

  React.useEffect(() => {
    const checkMobile = () => {
      setIsMobile(window.innerWidth <= 768);
    };
    checkMobile();
    window.addEventListener("resize", checkMobile);
    return () => {
      window.removeEventListener("resize", checkMobile);
    };
  }, []);

  const scaleDimensions = () => {
    return isMobile ? [0.7, 0.9] : [1.05, 1];
  };

  const rotate = useTransform(scrollYProgress, [0, 1], [20, 0]);
  const scale = useTransform(scrollYProgress, [0, 1], scaleDimensions());
  const translate = useTransform(scrollYProgress, [0, 1], [0, -100]);

  return (
    <div
      className={`h-[${height}rem] md:h-[${height}rem] flex items-center justify-center relative p-2 md:p-20 ` + className}
      ref={containerRef}
    >
      <div
        className={class2}
        style={{
          perspective: "1000px",
        }}
      >
        <Header translate={translate} titleComponent={titleComponent} />
        <div className="flex justify-around items-center">
          {orientation === "left" &&
            <><Card rotate={rotate} translate={translate} scale={scale}>
              {children}
            </Card><div className="flex justify-center items-center">
                <a href="#" className="block max-w-md text-justify p-6 pt-0 text-white shadow ">
                  <h5 className="mb-2 text-3xl font-bold tracking-tight text-white dark:text-white">{textHeader}</h5>
                  <p className="font-normal text-lg text-gray-200 dark:text-gray-400">{textDescription}</p>
                </a>
              </div></>
          }
          {orientation === "right" &&
            <>
              <div className="flex justify-center items-center">
                <a href="#" className="block max-w-md text-justify p-6 pt-0 text-white shadow ">
                  <h5 className="mb-2 text-3xl font-bold tracking-tight text-white dark:text-white">{textHeader}</h5>
                  <p className="font-normal text-lg text-gray-200 dark:text-gray-400">{textDescription}</p>
                </a>
              </div>
              <Card rotate={rotate} translate={translate} scale={scale}>
                {children}
              </Card>
            </>
          }


        </div>
      </div>
    </div>
  );
};

export const Header = ({ translate, titleComponent }: any) => {
  return (
    <motion.div
      style={{
        translateY: translate,
      }}
      className="div max-w-4xl mx-auto text-center"
    >
      {titleComponent}
    </motion.div>
  );
};

export const Card = ({
  rotate,
  scale,
  children,
}: {
  rotate: MotionValue<number>;
  scale: MotionValue<number>;
  translate: MotionValue<number>;
  children: React.ReactNode;
}) => {
  return (

    <motion.div
      style={{
        rotateX: rotate,
        scale,
        boxShadow:
          "0 0 #0000004d, 0 9px 20px #0000004a, 0 37px 37px #00000042, 0 84px 50px #00000026, 0 149px 60px #0000000a, 0 233px 65px #00000003",
      }}
      className="max-w-5xl justify-end items-center -mt-6 w-max h-max  border-4 border-[#6C6C6C] p-2 md:p-6 bg-[#222222] rounded-[30px] shadow-2xl"
    >
      <div className=" h-full w-full  overflow-hidden rounded-2xl bg-red-100 dark:bg-zinc-900 md:rounded-2xl md:p-4 ">
        {children}
      </div>
    </motion.div>
  );
};
