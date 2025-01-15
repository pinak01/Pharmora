import Footer from "@/components/footer";
import { Bento } from "@/components/grid";
import Landing from "@/components/landing";
import { HeroScrollDemo } from "@/components/scroll";
import { BackgroundBeams } from "@/components/ui/background-beams";
export default function Home() {
  return (
    <div className="w-full bg-[#000] h-full">
      <Landing/>
      <BackgroundBeams />
      <HeroScrollDemo/>
      <Bento/>
      <Footer></Footer>
    </div>
  );
}
