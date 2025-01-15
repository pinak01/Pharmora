import type { Metadata } from "next";
import "@/app/assets/globals.css"
import localFont from "next/font/local";
const akeila = localFont({ src: "../fonts/Akeila.otf" });

export const metadata: Metadata = {
  title: "Pharmora",
  description: "A comprehensive drug research platorm",
};

export default function RootLayout({
  children,
}: Readonly<{
  children: React.ReactNode;
}>) {
  return (
    <html lang="en">
      <body>
        {children}
      </body>
    </html>
  );
}
